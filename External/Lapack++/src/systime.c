
/*
*       double total =  systime(double &user, double &sys)
*
*   Returns the elapsed user and system time in seconds. Function
*   returns the elapsed time (> or = to  user+sys time).
*
*/

#ifndef WIN32
#include<sys/time.h>
#include<sys/resource.h>

#define tim(tv) ((double) (tv.tv_sec + 1.e-6 * tv.tv_usec))


double systime(user,system)
double *user, *system;
{
    static struct rusage r;
    static struct timeval tv;
    static struct timezone tz;
    
    getrusage(RUSAGE_SELF, &r);
    *user = tim(r.ru_utime);
    *system = tim(r.ru_stime);
    gettimeofday(&tv,&tz);
    return tim(tv);
}


#ifdef TEST_SYSTIME
#include <stdio.h>
main()
{
    double systime();
    double user[2], system[2], total[2];
    char s[80];
    double result= 0.0;
    int M;

for (;;)
  {
    printf("Enter # of loops to perform (MxM):"); fflush(stdout);
    scanf("%d", &M);
    total[0] = systime(&user[0], &system[0]);

    /* some useless computation */
    { int i; int j;
        for (i=0; i<M; i++)
            for (j=0; j<M; j++)
                result = result * 0.1 + ((double) i)*j /1e-6;
    }               

    total[1] = systime(&user[1], &system[1]);
    
    printf("Total: %8.2f,  User: %8.2f  System: %8.2f\n", 
        total[1]-total[0], user[1]-user[0], system[1]-system[0]);
  }
}
#endif //TEST_SYSTIME



#endif
