#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

enum MPI_tags { S_top, S_bottom, S_left, S_right};

double F (double x, double y)
{
    return 2*(x*x+ y*y)*(1 - 2*x*x*y*y)*exp(1 - x*x*y*y);
}

double phi (double x, double y)
{
    return exp(1 - x*x*y*y);
}
//структура для хранения данных о сетке
struct grid
{
	int rank, N1, N2, px, py, x_index_from, x_index_to, y_index_from, y_index_to;
	double *x_grid, *y_grid;
	double eps;
    bool top, bottom, left, right;

    double *send_message_top, *send_message_bottom, *send_message_left, *send_message_right;
    double *recv_message_top, *recv_message_bottom, *recv_message_left, *recv_message_right;
    MPI_Request* send_requests;
    MPI_Request* recv_requests;
    MPI_Comm comm;

	grid (int rank, MPI_Comm comm, double* x_grid, double* y_grid, int N1, int N2, int px, int py, double eps):
		rank (rank), comm (comm), x_grid (x_grid), y_grid (y_grid),
		send_message_top (NULL), send_message_bottom (NULL), send_message_left (NULL), send_message_right (NULL),
		recv_message_top (NULL), recv_message_bottom (NULL), recv_message_left (NULL), recv_message_right (NULL),
		send_requests (NULL), recv_requests (NULL),
		N1 (N1), N2 (N2),px (px), py (py), eps (eps),
		x_index_from (0), x_index_to (0), y_index_from (0), y_index_to (0),
		top (false), bottom (false), left (false), right (false) {
			int step1, step2;
			step1 = int(floor(1.0 * N1 / px));
			step2 = int(floor(1.0 * N2 / py));
			x_index_from = int(floor(1.0 * step1 * floor(1.0 * rank / py)));
			y_index_from = int(floor((double(rank % py)) * step2));

			if ((rank + 1) % py == 0)
				y_index_to = N2;
			else
				y_index_to = y_index_from + step2;

			if (rank >= (px-1)*py)
				x_index_to = N1;
			else
				x_index_to = x_index_from + step1;

			if (x_index_from == 0)
				top = true;
			if (y_index_from == 0)
				left = true;
			if (y_index_to == N1)
				right = true;
			if (x_index_to == N1)
				bottom = true;
		}

	int x_size()
	{
		if (bottom)
			return x_index_to - x_index_from + 1;
		else
			return x_index_to - x_index_from;
	}

	int y_size()
	{
		if (right)
			return y_index_to - y_index_from + 1;
		else
			return y_index_to - y_index_from;
	}

	int realij(int i, int j, int& grid_i, int& grid_j)
	{
		grid_i = x_index_from + i;
		grid_j = y_index_from + j;
	}

	double xval(int grid_i) {return x_grid[grid_i];}
	double yval(int grid_j) {return y_grid[grid_j];}
	double hx  (int grid_i) {return x_grid[grid_i+1] - x_grid[grid_i];}
	double hy  (int grid_j) {return y_grid[grid_j+1] - y_grid[grid_j];}
	int get_top_rank()      {return rank - py;}
	int get_bottom_rank()   {return rank + py;}
	int get_left_rank()     {return rank - 1;}
	int get_right_rank()    {return rank + 1;}
	bool is_border(int grid_i, int grid_j) {return ((grid_i == 0) || (grid_j == 0) || (grid_i == N1) || (grid_j == N2));}
};

void zone_sizes (int size, int& px, int& py);
//скалярное произведение
double braket (grid gr,  double* v1,  double* v2);
//вычисление дельты
void cdelta (grid gr,  double *func, double *delta_func, double f_top, double f_bottom, double f_left, double f_right, int i, int j, int grid_i, int grid_j);
//вычисление дельты для блока
void cadelta (grid gr, double* delta_func,  double* func);
//вычисление r
void cr (grid gr, double *r,  double *delta_p);
//вычисление g
void cg (grid gr, double *g, double *r, double alpha);
//вычисление p
void cp (grid gr, double *p, double* p_prev, double *g, double tau);
//вычисление нормы
double norm (grid gr, double *v1, double *v2);
//задание начального значения для p
void init_p_prev (grid gr, double* p_prev, double* phi);

//main
int main (int argc, char** argv)
{
    //переменные для подсчета времени
	double t1, t2;
    //граничные значения прямоугольника П = [A1,A2]X[B1,B2]
	double A1 = 0.0;
	double A2 = 2.0;
	double B1 = 0.0;
	double B2 = 2.0;
    //число узлов в сетке, передается в качестве параметров
	 int N1 = atoi (argv[1]);
	 int N2 = atoi (argv[2]);
	//погрешность
	double eps = 0.0001;
    //массивы, содержащие конкретные координаты для X и Y границы сетки соответственно
	double* x_grid = new double [N1+1];
	double* y_grid = new double [N2+1];
    //заполнение в соотетствии с формулой (8)
	for (int i = 0; i <= N1; i++) x_grid[i] = A2*i/N1 + A1*(1 - 1.0*i/N1);
	for (int j = 0; j <= N2; j++) y_grid[j] = B2*j/N2 + B1*(1 - 1.0*j/N2);
	//отвечают за стороны прямоугольника области
	int px, py;
    //id процесса и общее число процессов соответственно
	int rank, size;
	//стандартные MPI процедуры
	MPI_Init (&argc, &argv);
	//определяет id  (или rank) текущего процесса
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	//определение общего числа процессов
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	//запуск таймера
	t1 = MPI_Wtime();
    //определяет значения px, py по числу процессов
	zone_sizes(size, px, py);
    //отсеивание лишних процессов, на тот маловероятный случай, если их число не степень двойки
	if (rank < px*py)
	{
        //вывод информации о числе процессов и потоков
		if (rank == 0)
		{
	       	//#ifdef _OPENMP
	       	#ifdef _OPENMP
	        std::cout << "OpenMP Max-threads = " << omp_get_max_threads() << std::endl;
	        #endif
			std::cout << "px=" << px << " py=" << py << " size=" << size << std::endl;
	    }
        //инициализация участка сетки для каждого процесса
	    grid gr(rank, MPI_COMM_WORLD, x_grid, y_grid, N1, N2, px, py, eps);
        //массивы для основных переменных, участвующих в вычислениях
	    double* p =       new double [gr.x_size() * gr.y_size()];
	    double* p_prev =  new double [gr.x_size() * gr.y_size()];
	    double* g =       new double [gr.x_size() * gr.y_size()];
	    double* r =       new double [gr.x_size() * gr.y_size()];
	    double* delta_p = new double [gr.x_size() * gr.y_size()];
	    double* delta_r = new double [gr.x_size() * gr.y_size()];
	    double* delta_g = new double [gr.x_size() * gr.y_size()];

        //задание начальных значений
	    for (int i = 0; i < gr.x_size() * gr.y_size(); i++)
	    {
            g[i] = 0.0;
            r[i] = 0.0;
            delta_p[i] = 0.0;
            delta_g[i] = 0.0;
            delta_r[i] = 0.0;
        }
        //прочие промежуточные переменные, участвующие в вычислениях
	    double g_g = 1.0;
	    double dr_g = 1.0;
	    double r_g = 1.0;
	    double alpha = 0.0;
	    double tau = 0.0;
        //массив значений граничной функции фи на сетке
	    double* phi_on_grid = new double [gr.x_size() * gr.y_size()];
	    for (int i = 0; i < gr.x_size(); i++)
	    {
	    	for (int j = 0; j < gr.y_size(); j++)
	    	{
	    		int grid_i, grid_j;
	    		gr.realij(i, j, grid_i, grid_j);
	    		phi_on_grid[i*gr.y_size()+j] = phi(gr.xval(grid_i), gr.yval(grid_j));
	    	}
		}
		init_p_prev(gr, p_prev, phi_on_grid);
        //номер итерации
	    int n_iter = 1;
	    //основной цикл
	    while (true)
	    {
	    	cadelta(gr, delta_p, p_prev);
	    	cr(gr, r, delta_p);

	    	if (n_iter > 1)
	    	{
	    		cadelta(gr, delta_r, r);
	    		dr_g = braket(gr, delta_r, g);
	    		alpha = 1.0 * dr_g / g_g;
	    	}

	    	if (n_iter > 1)
	    		cg(gr, g, r, alpha);
	    	else
            	std::swap(g, r);

            cadelta(gr, delta_g, g);
            if (n_iter > 1)
            {
            	r_g = braket(gr, r, g);
            }
            else
            {
            	r_g = braket(gr, g, g);
            }

            g_g = braket(gr, delta_g, g);
	        tau = 1.0 * r_g / g_g;

	       	cp(gr, p, p_prev, g, tau);
	       	double norm_p = norm(gr, p, p_prev);
	       	if (rank == 0 && n_iter%100 == 1)
	       		printf("# iteration %d: norm: %f \n", n_iter, norm_p);
	       	if (norm_p < gr.eps && n_iter > 1)
            	break;

            swap(p, p_prev);
	    	n_iter += 1;
	    }

	}
    if (rank == 0)
    {
        //вывод затраченного времени
		t2 = MPI_Wtime();
        printf ("Elapsed time is %f\n", t2 - t1);
	}

    //завершение работы
	MPI_Finalize();
	return 0;
}


void zone_sizes (int size, int& px, int& py)
{
    if (size >= 512)
    {
        px = 16;
        py = 32;
    }
    else if (size >= 256)
    {
        px = 16;
        py = 16;
    }
    else if (size >= 128)
    {
        px = 8;
        py = 16;
    }
    else if (size >= 64)
    {
        px = 8;
        py = 8;
    }
    else if (size >= 32)
    {
        px = 4;
        py = 8;
    }
    else if (size >= 16)
    {
        px = 4;
        py = 4;
    }
    else if (size >= 8)
    {
        px = 2;
        py = 4;
    }
    else if (size >= 4)
    {
        px = 2;
        py = 2;
    }
    else if (size >= 2)
    {
        px = 1;
        py = 2;
    }
    else if (size >= 1)
    {
        px = 1;
        py = 1;
    }
}

double braket (grid gr,  double* v1,  double* v2)
{
	double res = 0.0;
	#pragma omp parallel for reduction(+:res)
	for (int i=0; i<gr.x_size(); i++)
	{
        for (int j=0; j<gr.y_size(); j++)
        {
        	int grid_i, grid_j;
	    	gr.realij(i, j, grid_i, grid_j);
        	if (not gr.is_border(grid_i, grid_j))
        	{
	        	double av_hx = (gr.hx(grid_i) + gr.hx(grid_i-1)) / 2.0;
	        	double av_hy = (gr.hy(grid_j) + gr.hy(grid_j-1)) / 2.0;
	            res += av_hx * av_hy * v1[i*gr.y_size()+j] * v2[i*gr.y_size()+j];
	        }
        }
    }

    double global_res = 0.0;
    MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_SUM, gr.comm);
    return global_res;
}

void cdelta (grid gr,  double *func, double *delta_func, double f_top, double f_bottom, double f_left, double f_right, int i, int j, int grid_i, int grid_j)
{
	double h_i_1 = gr.hx(grid_i-1);
	double h_i = gr.hx(grid_i);
	double h_j_1 = gr.hy(grid_j-1);
	double h_j = gr.hy(grid_j);
	double av_hx = (h_i + h_i_1) / 2.0;
	double av_hy = (h_j + h_j_1) / 2.0;
	double f_curr = func[i*gr.y_size()+j];
	delta_func[i*gr.y_size()+j] =
		(1.0 / av_hx) * ((f_curr - f_top) / h_i_1 - (f_bottom - f_curr) / h_i) +
		(1.0 / av_hy) * ((f_curr - f_left) / h_j_1 - (f_right - f_curr) / h_j);
}

void cadelta (grid gr, double* delta_func,  double* func)
{
	//вычисление внутренних точек
	int i, j;
	#pragma omp parallel for
	for (i = 1; i < gr.x_size()-1; i++)
	{
    	for (j = 1; j < gr.y_size()-1; j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		cdelta(gr, func, delta_func, func[(i-1)*gr.y_size()+j], func[(i+1)*gr.y_size()+j],
                            func[i*gr.y_size()+j-1], func[i*gr.y_size()+j+1], i, j, grid_i, grid_j);
    	}
	}

	if (gr.send_message_top == NULL)
		gr.send_message_top = new double [gr.y_size()];
	if (gr.send_message_bottom == NULL)
		gr.send_message_bottom = new double [gr.y_size()];
	if (gr.send_message_left == NULL)
		gr.send_message_left = new double [gr.x_size()];
	if (gr.send_message_right == NULL)
		gr.send_message_right = new double [gr.x_size()];

	if (gr.recv_message_top == NULL)
		gr.recv_message_top = new double [gr.y_size()];
	if (gr.recv_message_bottom == NULL)
		gr.recv_message_bottom = new double [gr.y_size()];
	if (gr.recv_message_left == NULL)
		gr.recv_message_left = new double [gr.x_size()];
	if (gr.recv_message_right == NULL)
		gr.recv_message_right = new double [gr.x_size()];

	if (gr.send_requests == NULL)
		gr.send_requests = new MPI_Request [4];
	if (gr.recv_requests == NULL)
		gr.recv_requests = new MPI_Request [4];

    //граничные значения области пересылаются соотв. процессам (подготовка)
	for (int j = 0; j < gr.y_size(); j++)
	{
		gr.send_message_top[j] = func[0*gr.y_size()+j];
		gr.send_message_bottom[j] = func[(gr.x_size()-1)*gr.y_size()+j];
    }
	for (int i = 0; i < gr.x_size(); i++)
	{
		gr.send_message_left[i] = func[i*gr.y_size()+0];
		gr.send_message_right[i] = func[i*gr.y_size()+gr.y_size()-1];
    }

	int status;
	int send_count = 0;
	if (not gr.top)
	{
		MPI_Isend(gr.send_message_top, gr.y_size(), MPI_DOUBLE,
			gr.get_top_rank(), S_top, gr.comm, &(gr.send_requests[send_count]));
		send_count++;
	}
	if (not gr.bottom)
	{
		MPI_Isend(gr.send_message_bottom, gr.y_size(), MPI_DOUBLE,
			gr.get_bottom_rank(), S_bottom, gr.comm, &(gr.send_requests[send_count]));
		send_count++;
	}
	if (not gr.left)
	{
		MPI_Isend(gr.send_message_left, gr.x_size(), MPI_DOUBLE,
			gr.get_left_rank(), S_left, gr.comm, &(gr.send_requests[send_count]));
		send_count++;
	}
	if (not gr.right)
	{
		MPI_Isend(gr.send_message_right, gr.x_size(), MPI_DOUBLE,
			gr.get_right_rank(), S_right, gr.comm, &(gr.send_requests[send_count]));
		send_count++;
	}

	int recv_count = 0;
	if (not gr.top)
	{
		MPI_Irecv(gr.recv_message_top, gr.y_size(), MPI_DOUBLE,
			gr.get_top_rank(), S_bottom, gr.comm, &(gr.recv_requests[recv_count]));
		recv_count++;
	}
	if (not gr.bottom)
	{
		MPI_Irecv(gr.recv_message_bottom, gr.y_size(), MPI_DOUBLE,
			gr.get_bottom_rank(), S_top, gr.comm, &(gr.recv_requests[recv_count]));
		recv_count++;
	}
	if (not gr.left) {
		MPI_Irecv(gr.recv_message_left, gr.x_size(), MPI_DOUBLE,
			gr.get_left_rank(), S_right, gr.comm, &(gr.recv_requests[recv_count]));
		recv_count++;
	}
	if (not gr.right)
	{
		MPI_Irecv(gr.recv_message_right, gr.x_size(), MPI_DOUBLE,
			gr.get_right_rank(), S_left, gr.comm, &(gr.recv_requests[recv_count]));
		recv_count++;
	}

	MPI_Waitall(recv_count, gr.recv_requests, MPI_STATUS_IGNORE);
    MPI_Waitall(send_count, gr.send_requests, MPI_STATUS_IGNORE);

    //вычисление краев
    if (not gr.top)
    {
    	int i = 0;
    	for (int j = 1; j < gr.y_size()-1; j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		cdelta(gr, func, delta_func, gr.recv_message_top[j], func[(i+1)*gr.y_size()+j],
                            func[i*gr.y_size()+j-1], func[i*gr.y_size()+j+1], i, j, grid_i, grid_j);
    	}
    }

	if (not gr.bottom)
	{
    	int i = gr.x_size()-1;
    	for (int j = 1; j < gr.y_size()-1; j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		cdelta(gr, func, delta_func, func[(i-1)*gr.y_size()+j], gr.recv_message_bottom[j],
                            func[i*gr.y_size()+j-1], func[i*gr.y_size()+j+1], i, j, grid_i, grid_j);
    	}
    }

    if (not gr.left)
    {
    	int j = 0;
    	for (int i = 1; i < gr.x_size()-1; i++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		cdelta(gr, func, delta_func, func[(i-1)*gr.y_size()+j], func[(i+1)*gr.y_size()+j],
                            gr.recv_message_left[i], func[i*gr.y_size()+j+1], i, j, grid_i, grid_j);
    	}
    }

    if (not gr.right)
    {
    	int j = gr.y_size()-1;
    	for (int i = 1; i < gr.x_size()-1; i++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		cdelta(gr, func, delta_func, func[(i-1)*gr.y_size()+j], func[(i+1)*gr.y_size()+j],
                            func[i*gr.y_size()+j-1], gr.recv_message_right[i], i, j, grid_i, grid_j);
    	}
    }

    //вычисление угловых точек
	i = 0; j = 0;
	if (not gr.top && not gr.left)
	{
		int grid_i, grid_j;
    	gr.realij(i, j, grid_i, grid_j);
    	cdelta(gr, func, delta_func, gr.recv_message_top[j], func[(i+1)*gr.y_size()+j],
                        gr.recv_message_left[i], func[i*gr.y_size()+j+1], i, j, grid_i, grid_j);
	}

	i = 0; j = gr.y_size()-1;
	if (not gr.top && not gr.right)
	{
		int grid_i, grid_j;
    	gr.realij(i, j, grid_i, grid_j);
    	cdelta(gr, func, delta_func, gr.recv_message_top[j], func[(i+1)*gr.y_size()+j],
                        func[i*gr.y_size()+j-1], gr.recv_message_right[i], i, j, grid_i, grid_j);
	}

	i = gr.x_size()-1; j = 0;
	if (not gr.bottom && not gr.left)
	{
		int grid_i, grid_j;
    	gr.realij(i, j, grid_i, grid_j);
    	cdelta(gr, func, delta_func, func[(i-1)*gr.y_size()+j], gr.recv_message_bottom[j],
                        gr.recv_message_left[i], func[i*gr.y_size()+j+1], i, j, grid_i, grid_j);
	}

	i = gr.x_size()-1; j = gr.y_size()-1;
	if (not gr.bottom && not gr.right)
	{
		int grid_i, grid_j;
    	gr.realij(i, j, grid_i, grid_j);
    	cdelta(gr, func, delta_func, func[(i-1)*gr.y_size()+j], gr.recv_message_bottom[j],
                        func[i*gr.y_size()+j-1], gr.recv_message_right[i], i, j, grid_i, grid_j);
	}
}

void cr (grid gr, double *r,  double *delta_p)
{
	int i, j;
	#pragma omp parallel for
	for (i = 0; i < gr.x_size(); i++)
	{
    	for (j = 0; j < gr.y_size(); j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		if (gr.is_border(grid_i, grid_j))
				r[i*gr.y_size()+j] = 0.0;
    		else
    			r[i*gr.y_size()+j] = delta_p[i*gr.y_size()+j] - F(gr.xval(grid_i), gr.yval(grid_j));
    	}
	}
}

void cg (grid gr, double *g, double *r, double alpha)
{
	int i, j;
	#pragma omp parallel for
	for (i = 0; i < gr.x_size(); i++)
	{
    	for (j = 0; j < gr.y_size(); j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		g[i*gr.y_size()+j] = r[i*gr.y_size()+j] - alpha * g[i*gr.y_size()+j];
    	}
	}
}

void cp (grid gr, double *p, double* p_prev, double *g, double tau)
{
	int i, j;
	#pragma omp parallel for
	for (i = 0; i < gr.x_size(); i++)
	{
    	for (j = 0; j < gr.y_size(); j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		p[i*gr.y_size()+j] = p_prev[i*gr.y_size()+j] - tau * g[i*gr.y_size()+j];
    	}
	}
}

double norm (grid gr, double *v1, double *v2)
{
	double norm_shared = 0.0;
    #pragma omp parallel for reduction(+:norm_shared)
	for (int i = 0; i < gr.x_size(); i++)
	{
        for (int j = 0; j < gr.y_size(); j++)
        {
        	int grid_i, grid_j;
	    	gr.realij(i, j, grid_i, grid_j);
        	if (not gr.is_border(grid_i, grid_j))
        	{
	        	double av_hx = (gr.hx(grid_i) + gr.hx(grid_i-1)) / 2.0;
	        	double av_hy = (gr.hy(grid_j) + gr.hy(grid_j-1)) / 2.0;
                double tmp = (v1[i*gr.y_size()+j] - v2[i*gr.y_size()+j]);
	            norm_shared += av_hx * av_hy * tmp*tmp;
	        }
        }
    }
	double global_norm = 0.0;
	MPI_Allreduce(&norm_shared, &global_norm, 1, MPI_DOUBLE, MPI_SUM, gr.comm);
    return sqrt(global_norm);
}

void init_p_prev (grid gr, double* p_prev, double* phi)
{
	int i, j;
	#pragma omp parallel for
	for (i = 0; i < gr.x_size(); i++)
	{
    	for (j = 0; j < gr.y_size(); j++)
    	{
    		int grid_i, grid_j;
    		gr.realij(i, j, grid_i, grid_j);
    		if (not gr.is_border(grid_i, grid_j))
    		{
                p_prev[i*gr.y_size()+j] = 0.0;
            }
            else
            {
                p_prev[i*gr.y_size()+j] = phi[i*gr.y_size()+j];
            }
		}
	}
}


//end

