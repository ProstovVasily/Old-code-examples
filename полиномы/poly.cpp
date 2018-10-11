#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <mutex>
#include <thread>
#include <condition_variable>



using std::ifstream;
using std::ofstream;
using std::string;
using std::to_string;
using std::vector;
using std::queue;
using std::mutex;
using std::thread;
using std::unique_lock;
using std::condition_variable;
using std::ref;
using std::cout;
using std::endl;

//add here some mini class to store polynomes.
//parsing string into massive of monomes will be done inside the constructor
class polynoms
{
//who cares, lets make everything public for simplicity! (just kidding, i should delete it later)
public:
    //actually parsing string version to our representation,
    polynoms (string s, int n);
    //not important
    polynoms () {nv = 0;}
    // calculates polynom's value given that x is a vector of 1 and 0 with nv length
    int P(vector<int> x);
    //temporal, for checking polynomps
    void print();

private:
    //number of variables
    int nv;
    // whether 1 in polynom or not
    bool one;
    //vector of monomes, each holding information of whether x_i in it or not
    vector <vector <int>> poly;
};

//structure to be sent to queue
struct wrap
{
public:
    //mutex, responsible for information about wheter this part already have been processed or not
    mutex& m;
    //will contain results of calculations
    string& res;
    //polynom used in calcucations
    polynoms poly;
    //number of variavles in polynom
    int nv;
    //shift from the begining
    int shift;

    wrap (mutex& m_, string& res_, polynoms poly_, int nv_, int shift_): m(m_), res(res_), poly(poly_), nv(nv_), shift(shift_){};

};


//queue with some parts to process
class blocking_queue
{
public:
    //take front element from queue
    wrap pop();
    //push element into queue
    void push(const wrap& item);

private:
    //actual queue
    queue<wrap> q;
    //used for the situation with empty queue
    mutex m;
    condition_variable con;
};


//reading & writing from file input_i to output_i. nv represents number of variables and ns represents number of solvers.
void reader (int i, int nv, int ns, blocking_queue& q);

//solving some stuff, calculating part of the polynom's vector.
void solver (blocking_queue& q);

//convert number x < 2^n into vector of 1 and 0
vector<int> int_to_vector (int x, int n);

int main (int argc, char** argv)
{
    int ns = 4, nr = 1, nv=3;
    if (argc > 1) nv = atoi(argv[1]);
    if (argc > 2) nr = atoi(argv[2]);
    if (argc > 3) ns = atoi(argv[3]);
    ns = ns > 2 ? ns : 2;
    blocking_queue q;
    thread readers[nr];
    thread solvers[ns];

    for (int i = 0; i < nr; i++) readers[i] = thread(reader, i+1, nv, ns, ref(q));
    for (int i = 0; i < ns; i++) solvers[i] = thread(solver, ref(q));
    for (int i = 0; i < nr; i++) readers[i].join();
    for (int i = 0; i < ns; i++) solvers[i].detach();

    return 0;
}

void reader (int i, int nv, int ns, blocking_queue& q)
{
    ifstream inp ("input_"  + to_string(i) + ".txt");
    ofstream out ("output_" + to_string(i) + ".txt");
    string s;
    mutex solv[ns];
    string res[ns];
    int len = (2 << (nv-1))/(ns-1);
    for (int i = 0; i < ns-1; i++) res[i] = string(len, '0');
    res[ns-1] = string ((2 << (nv-1)) - len*(ns-1), '0');
    while (inp >> s)//change to read lines rather than strings until ' '//later
    {
        polynoms pol(s, nv);
        string ans= "(";
        for (int i = 0; i < ns; i++)
        {
            solv[i].lock();
            q.push(wrap(solv[i], res[i], pol, nv, len*i));
        }
        for (int i = 0; i < ns; i++)
        {
            solv[i].lock();
            solv[i].unlock();
            ans += res[i];
        }
        ans += ")";
        out << ans << std::endl;
    }
}

polynoms::polynoms(string s, int n)
{
    nv = n;
    one = false;
    //here string is parsed into polynom;
    char c;
    vector <int> monom (nv, 0);
    int tmp = 0;
    bool flag = false;
    for (int i = 0; i < s.length(); i++)
    {
        //
        c = s[i];
        if (c == 'x' && !flag)
        {
            flag = true;
            continue;
        }
        if (c == 'x' && flag)
        {
            if (tmp == 0) continue;
            if (tmp > nv)
            {
                tmp = 0;
                continue;
            }
            tmp -= 1;
            monom[tmp] = 1;
            tmp = 0;
            continue;
        }
        if (c == '1' && !flag)
        {
            //change inner flag for 1
            one = !one;
            continue;
        }
        if (c >= '0' && c <= '9' && flag)
        {
            tmp = tmp*10 + (c - '0');
            continue;
        }
        if (c == '+' && flag)
        {
            //
            if (tmp > nv)
            {
                tmp = 0;
                poly.push_back(monom);
                monom = vector<int> (nv, 0);
                flag = false;
                continue;
            }
            tmp -= 1;
            monom[tmp] = 1;
            tmp = 0;
            //append
            poly.push_back(monom);
            monom = vector<int> (nv, 0);
            flag = false;
        }
    }
    if (flag)
    {
        if (tmp <= nv)
        {
            tmp -= 1;
            monom[tmp] = 1;
            tmp = 0;
        }
        poly.push_back(monom);
    }

}

void polynoms::print ()
{
    for (int i = 0; i < poly.size(); i++)
    {
        for (int j = 0; j < nv; j++)
        {
            std::cout << poly[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

int polynoms::P(vector<int> x)
{
    int ans = 0;
    if (one) ans += 1;
    int tmp = 0;
    for (int i = 0; i < poly.size(); i++)
    {
        tmp = 0;
        for (int j = 0; j < nv; j++)
        {
            if (poly[i][j] == 0) continue;
            if (x[j] == 1) tmp = 1;
            else
            {
                tmp = 0;
                break;
            }
        }
        ans = (ans + tmp)%2;
    }
    return ans;
}

vector<int> int_to_vector(int x, int n)
{
    //if x >  2^n then whatever
    vector<int> ans;

    for (int i = 0; i < n; i++)
    {
        ans.push_back(x%2);
        x = x/2;
    }
    std::reverse(ans.begin(), ans.end());
    return ans;
}

void solver (blocking_queue& q)
{
    int tmp = 0;
    while (true)
    {
        wrap w = q.pop();
        for (int i = w.shift; i < w.res.length() + w.shift; i++)
        {
            tmp = w.poly.P(int_to_vector(i, w.nv));
            w.res[i - w.shift] = '0' + tmp;
        }
        w.m.unlock();
    }
}

wrap blocking_queue::pop()
{
    unique_lock<mutex> mlock(m);
    while (q.empty())
    {
      con.wait(mlock);
    }
    wrap res = q.front();
    q.pop();
    return res;
}


void blocking_queue::push(const wrap& x)
{
    unique_lock<mutex> mlock(m);
    q.push(x);
    mlock.unlock();
    con.notify_one();
}








