#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>
#define DEFAULT_COMM_FOLDER "/bin/"

char ** read_com (char* last);
int spec (char tmp);
void free_com (char** com);
int run_com (char** com, int red_inp, int red_out);
int bg_com (char** com);

int main ()
{
    int sz_i = 10;
    char res = 0;
    int is_conv = 0;
    int k;
    int add;
    char path [20];
    //int fd [10] [2];
    int f;
    char** tmp;
    int (*fd)[2] = calloc (sz_i, sizeof (int[2]));
    char*** test = calloc (sz_i, sizeof (char**));
    int i = 0;
    for (i = 0;;i++)
    {
        if (i >= sz_i - 1)
        {
            sz_i *= 2;
            test = realloc (test, sz_i*sizeof (char**));
            fd = realloc (fd, sz_i*sizeof (int[2]));
        }
        test[i] = read_com (&res);
        if (res == '<')
        {
            scanf ("%s", path);
            k = open (path, O_RDONLY);
            f = run_com (test[i], k, 0);
            waitpid (f, 0, 0);
        }
        else
        if (is_conv)
        {
            if (res != '|')
            {
                k = 0;
                if (res == '>')
                {
                    res = getchar ();
                    if (res != '>')
                    {
                        ungetc (res, stdin);
                        add = O_TRUNC;
                    }
                    else add = O_APPEND;
                    scanf ("%s", path);
                    k = open (path, O_WRONLY | add | O_CREAT);
                }
                is_conv = 0;
                f = run_com (test[i], fd[i-1][0], k);
                close (fd[i-1][0]);
                waitpid (f, 0, 0);
                if ((res == '!') || (res == ')')) break;
                if (res == '(') main ();
            }
            else
            {
                res = getchar ();
                if (res == '|')
                {
                    pipe (fd[i]);
                    run_com (test[i], fd[i-1][0], fd[i][1]);
                    close (fd[i][1]);
                    close (fd[i-1][0]);
                    waitpid (f, &k, 0);
                    if (k == 0)
                    {
                        tmp = read_com (&res);
                        free_com (tmp);
                    }
                }
                else
                {
                    ungetc (res, stdin);
                    pipe (fd[i]);
                    run_com (test[i], fd[i-1][0], fd[i][1]);
                    close (fd[i][1]);
                    close (fd[i-1][0]);
                }
            }
        }
        else
        if (res == '&')
        {
            res = getchar ();
            if (res != '&')
            {
                ungetc (res, stdin);
                bg_com (test[i]);
            }
            else
            {
                //if &&
                //not done yet
                f = run_com (test[i], 0, 0);
                waitpid (f, &k, 0);
                if (k != 0)
                {
                    tmp = read_com (&res);
                    free_com (tmp);
                }
            }
        }
        else
        if (res == '>')
        {
            res = getchar ();
            if (res != '>')
            {
                ungetc (res, stdin);
                add = O_TRUNC;
            }
            else add = O_APPEND;
            scanf ("%s", path);
            k = open (path, O_WRONLY | add | O_CREAT, 0777);
            f = run_com (test[i], 0, k);
            waitpid (f, 0, 0);
        }
        else
        if (res == ';')
        {
            f = run_com (test[i], 0, 0);
            waitpid (f, 0, 0);
        }
        else
        if (res == '|')
        {
            res = getchar ();
            if (res != '|')
            {
                ungetc (res, stdin);
                pipe (fd[i]);
                run_com (test[i], 0, fd[i][1]);
                is_conv = 1;
                close (fd[i][1]);
            }
            else
            {
                //not done yet
                f = run_com (test[i], 0, 0);
                waitpid (f, &k, 0);
                if (k == 0)
                {
                    tmp = read_com (&res);
                    free_com (tmp);
                }
            }
        }
        else
        if (res == '\n')
        {
            f = run_com (test[i], 0, 0);
            waitpid (f, 0, 0);
        }
        else
        if (res == '!')
        {
            f = run_com (test[i], 0, 0);
            waitpid (f, 0, 0);
            break;
        }
        if (res == ')')
        {
            run_com (test[i], 0, 0);
            break;
        }
        if (res == '(')
        {
            run_com (test[i], 0, 0);
            main ();
        }
        //if (!is_conv) wait (0);
    }
    printf ("111\n");
    for (k = 0; k < i; k++)
    {
        free_com (test[k]);
        wait (0);
    }
    free (test);
    return 0;
}

char** read_com (char* last)
{
    char ** res;
    int i = 0;
    int j = 0;
    int sz_i = 10;
    int sz_j = 10;
    char tmp;
    tmp = getchar ();
    res = calloc (sz_i, sizeof (char*));
    while (tmp == ' ') tmp = getchar();
    if (!spec(tmp))
    {
        *last = tmp;
        free (res);
        return NULL;
    }
    while (spec(tmp))
    {
        while (tmp == ' ') tmp = getchar();
        if (i >= sz_i - 1)
        {
            sz_i *= 2;
            res = realloc (res, sz_i*sizeof (char*));
        }
        if (!spec(tmp))
        {
            res[i] = NULL;
            break;
        }
        j = 0;
        res[i] = calloc (sz_i, sizeof (char));
        while (spec(tmp) && (tmp != ' '))
        {
            if (j >= sz_j - 1)
            {
                sz_i *= 2;
                res[i] = realloc (res[i], sz_j*sizeof (char));
            }
            res[i][j] = tmp;
            tmp = getchar ();
            j++;
        }
        res[i][j] = '\0';
        i++;
        if (!spec(tmp))
        {
            res[i] = NULL;
            break;
        }
    }
    *last = tmp;
    return res;
}

int spec (char temp)
{
    if (temp == '&') return 0;
    if (temp == '|') return 0;
    if (temp == ';') return 0;
    if (temp == '>') return 0;
    if (temp == '<') return 0;
    if (temp == '(') return 0;
    if (temp == ')') return 0;
    if (temp == '!') return 0;
    if (temp == '\n') return 0;
    return 1;
}

void free_com (char** com)
{
    int i = 0;
    if (com == NULL) return;
    while (com[i] != NULL) free (com[i++]);
    free (com);
}

int run_com (char** com, int red_inp, int red_out)
{
    unsigned char isFile =0x8;
    DIR* dir;
    int f;
    struct dirent* dir_info;
    char* newpath = NULL;
    if (com == NULL) return 0;
    dir = opendir (DEFAULT_COMM_FOLDER);
    while (dir_info = readdir (dir))
    {
        if (dir_info->d_type == isFile)
            if (strcmp (dir_info->d_name, com[0]) == 0)
            {
                newpath = calloc (strlen (DEFAULT_COMM_FOLDER) + strlen (dir_info->d_name) + 1, sizeof (char));
                strcpy (newpath, DEFAULT_COMM_FOLDER);
                strcat (newpath, dir_info->d_name);
                if (!(f = fork()))
                {
                    if (red_inp)
                    {
                        close (0);
                        dup (red_inp);
                    }
                    if (red_out)
                    {
                        close (1);
                        dup (red_out);
                    }
                    execvp (newpath, com);
                    exit (0);
                }
                return f;
            };
    }
    if (!(f = fork()))
    {
        if (red_inp)
        {
            close (0);
            dup (red_inp);
        }
        if (red_out)
        {
            close (1);
            dup (red_out);
        }
        if (execvp (com[0], com) == -1)
        {
            printf ("no such command or file '%s'\n", com[0]);
            exit (1);
        }
        exit (0);
    }
    return f;
}

int bg_com (char** com)
{
    //add copying test to smth temporary;
    if (!fork())
    {
        if (!fork())
        {
            signal (SIGINT, SIG_IGN);
            close (0);
            close (1);
            open ("/dev/null", O_RDONLY);
            open ("/dev/null", O_WRONLY);
            run_com (com, 0, 0);
            exit (0);
        }
        exit (0);
    }
    return 0;
}

