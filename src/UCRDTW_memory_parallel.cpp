/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected © 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define dist(x, y) ((x - y) * (x - y))

using namespace std;

/// Data structure for sorting the query
typedef struct Index
{
    float value;
    int index;
} Index;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
struct deque
{
    int *dq;
    int size, capacity;
    int f, r;
};

typedef struct transdata
{
    long int startnumber, stopnumber;
    float *u, *l, *qo, *uo, *lo, *tz, *cb, *cb1, *cb2, *u_d, *l_d;
    float *buffer, *u_buff, *l_buff;
    float *t, *q;
    float *rawfile;
    int *order;
    float *bsf;
    int m = -1, r = -1, lengthSeriesDataset = 0, numberSeriesDataset = 0;
    bool bNormalized;
} transdata;
/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void *b)
{
    Index *x = (Index *)a;
    Index *y = (Index *)b;
    return abs(y->value) - abs(x->value); // high to low
}

/// Initial the queue at the begining step of envelop calculation
void init(deque *d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *)malloc(sizeof(int) * d->capacity);
    d->f = 0;
    d->r = d->capacity - 1;
}

/// Destroy the queue
void destroy(deque *d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque *d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity - 1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque *d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity - 1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque *d)
{
    d->r = (d->r + 1) % d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque *d)
{
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity - 1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque *d)
{
    int aux = (d->r + 1) % d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque *d)
{
    return d->size == 0;
}

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(float *t, int len, int r, float *l, float *u)
{
    struct deque du, dl;

    init(&du, 2 * r + 2);
    init(&dl, 2 * r + 2);

    push_back(&du, 0);
    push_back(&dl, 0);

    for (int i = 1; i < len; i++)
    {
        if (i > r)
        {
            u[i - r - 1] = t[front(&du)];
            l[i - r - 1] = t[front(&dl)];
        }
        if (t[i] > t[i - 1])
        {
            pop_back(&du);
            while (!empty(&du) && t[i] > t[back(&du)])
                pop_back(&du);
        }
        else
        {
            pop_back(&dl);
            while (!empty(&dl) && t[i] < t[back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    for (int i = len; i < len + r + 1; i++)
    {
        u[i - r - 1] = t[front(&du)];
        l[i - r - 1] = t[front(&dl)];
        if (i - front(&du) >= 2 * r + 1)
            pop_front(&du);
        if (i - front(&dl) >= 2 * r + 1)
            pop_front(&dl);
    }
    destroy(&du);
    destroy(&dl);
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
float lb_kim_hierarchy(float *t, float *q, int j, int len, float mean, float std, _Bool bNormalized, float bsf = FLT_MAX)
{
    /// 1 point at front and back
    float d, lb;

    float x0;
    float y0;

    if (bNormalized)
    {
        x0 = (t[j] - mean) / std;
        y0 = (t[(len - 1 + j)] - mean) / std;
    }
    else
    {
        x0 = t[j];
        y0 = t[(len - 1 + j)];
    }

    lb = dist(x0, q[0]) + dist(y0, q[len - 1]);
    if (lb >= bsf)
        return lb;

    /// 2 points at front
    float x1;
    if (bNormalized)
        x1 = (t[(j + 1)] - mean) / std;
    else
        x1 = t[(j + 1)];

    d = min(dist(x1, q[0]), dist(x0, q[1]));
    d = min(d, dist(x1, q[1]));
    lb += d;
    if (lb >= bsf)
        return lb;

    /// 2 points at back

    float y1;
    if (bNormalized)
        y1 = (t[(len - 2 + j)] - mean) / std;
    else
        y1 = t[(len - 2 + j)];

    d = min(dist(y1, q[len - 1]), dist(y0, q[len - 2]));
    d = min(d, dist(y1, q[len - 2]));
    lb += d;
    if (lb >= bsf)
        return lb;

    /// 3 points at front
    float x2;
    if (bNormalized)
        x2 = (t[(j + 2)] - mean) / std;
    else
        x2 = t[(j + 2)];

    d = min(dist(x0, q[2]), dist(x1, q[2]));
    d = min(d, dist(x2, q[2]));
    d = min(d, dist(x2, q[1]));
    d = min(d, dist(x2, q[0]));
    lb += d;
    if (lb >= bsf)
        return lb;

    /// 3 points at back
    float y2;
    if (bNormalized)
        y2 = (t[(len - 3 + j)] - mean) / std;
    else
        y2 = t[(len - 3 + j)];

    d = min(dist(y0, q[len - 3]), dist(y1, q[len - 3]));
    d = min(d, dist(y2, q[len - 3]));
    d = min(d, dist(y2, q[len - 2]));
    d = min(d, dist(y2, q[len - 1]));
    lb += d;

    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
float lb_keogh_cumulative(int *order, float *t, float *uo, float *lo, float *cb, int j, int len, float mean, float std, _Bool bNormalized, float best_so_far = FLT_MAX)
{
    float lb = 0;
    float x, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        if (bNormalized)
            x = (t[(order[i] + j)] - mean) / std;
        d = 0;
        if (x > uo[i])
            d = dist(x, uo[i]);
        else if (x < lo[i])
            d = dist(x, lo[i]);
        else
            d = 0;
        if (t[(order[i] + j)] > uo[i])
            d = dist(t[(order[i] + j)], uo[i]);
        else if (t[(order[i] + j)] < lo[i])
            d = dist(t[(order[i] + j)], lo[i]);

        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
float lb_keogh_data_cumulative(int *order, float *tz, float *qo, float *cb, float *l, float *u, int len, float mean, float std, _Bool bNormalized, float best_so_far = FLT_MAX)
{
    float lb = 0;
    float uu, ll, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        d = 0;
        if (bNormalized)
        {
            uu = (u[order[i]] - mean) / std;
            ll = (l[order[i]] - mean) / std;

            if (qo[i] > uu)
                d = dist(qo[i], uu);
            else
            {
                if (qo[i] < ll)
                    d = dist(qo[i], ll);
            }
        }
        else
        {
            if (qo[i] > u[order[i]])
                d = dist(qo[i], u[order[i]]);
            else
            {
                if (qo[i] < l[order[i]])
                    d = dist(qo[i], l[order[i]]);
            }
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
float dtw(float *A, float *B, float *cb, int m, int r, int *pointChecked, float bsf = FLT_MAX)
{

    float *cost;
    float *cost_prev;
    float *cost_tmp;
    int i, j, k;
    float x, y, z, min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (float *)malloc(sizeof(float) * (2 * r + 1));
    for (k = 0; k < 2 * r + 1; k++)
        cost[k] = FLT_MAX;

    cost_prev = (float *)malloc(sizeof(float) * (2 * r + 1));
    for (k = 0; k < 2 * r + 1; k++)
        cost_prev[k] = FLT_MAX;

    for (i = 0; i < m; i++)
    {
        k = max(0, r - i);
        min_cost = FLT_MAX;

        for (j = max(0, i - r); j <= min(m - 1, i + r); j++, k++)
        {
            /// Initialize all row and column
            if ((i == 0) && (j == 0))
            {
                cost[k] = dist(A[0], B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0))
                y = FLT_MAX;
            else
                y = cost[k - 1];
            if ((i - 1 < 0) || (k + 1 > 2 * r))
                x = FLT_MAX;
            else
                x = cost_prev[k + 1];
            if ((i - 1 < 0) || (j - 1 < 0))
                z = FLT_MAX;
            else
                z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min(min(x, y), z) + dist(A[i], B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i + r < m - 1 && min_cost + cb[i + r + 1] >= bsf)
        {
            free(cost);
            free(cost_prev);
            return min_cost + cb[i + r + 1];
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
        (*pointChecked)++;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    float final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

/// Print function for debugging
void printArray(float *x, int len)
{
    for (int i = 0; i < len; i++)
        printf(" %6.2lf", x[i]);
    printf("\n");
}

/// If expected error happens, teminated the program.
void error(int id)
{
    if (id == 1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if (id == 2)
        printf("ERROR : File not Found!!!\n\n");
    else if (id == 3)
        printf("ERROR : Can't create Output File!!!\n\n");
    else if (id == 4)
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
        printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   [R (number of points)] [number series datasets] [length of series dataset] [ 1 normalized - 0 non-normalized] [queries number] [thread number] [optional: epsilon] \n\n");
        printf("For example for 1-NN  :  UCR_DTW.exe  data.bin   query.txt   128   5   4000    256  1 10 8 \n");
        printf("For example for epsilon-range (epsilon=5):  UCR_DTW.exe  data.bin   query.txt   128   5   4000  256  1 10 8 5 \n");
    }
    exit(1);
}
void *workerthread(void *transferdata)
{
    long int startnumber = ((transdata *)transferdata)->startnumber;
    long int stopnumber = ((transdata *)transferdata)->stopnumber;
    float *rawfile = ((transdata *)transferdata)->rawfile;
    float ex, ex2, mean, std;
    bool bNormalized;
    float bsf = FLT_MAX, dist = FLT_MAX, lb_kim = FLT_MAX, lb_k = FLT_MAX, lb_k2 = FLT_MAX;
    float *u, *l, *qo, *uo, *lo, *tz, *cb, *cb1, *cb2, *u_d, *l_d;
    float *buffer, *u_buff, *l_buff;
    float *t, *q;
    int m = -1, r = -1, lengthSeriesDataset = 256;
    int *order;
    float d;

    rawfile = ((transdata *)transferdata)->rawfile;
    q = ((transdata *)transferdata)->q;
    u = ((transdata *)transferdata)->u;
    l = ((transdata *)transferdata)->l;
    uo = ((transdata *)transferdata)->uo;
    lo = ((transdata *)transferdata)->lo;
    qo = ((transdata *)transferdata)->qo;
    r = ((transdata *)transferdata)->r;
    m = ((transdata *)transferdata)->m;
    bNormalized = ((transdata *)transferdata)->bNormalized;
    order = ((transdata *)transferdata)->order;
    cb = (float *)malloc(sizeof(float) * m);
    cb1 = (float *)malloc(sizeof(float) * m);
    cb2 = (float *)malloc(sizeof(float) * m);
    u_d = (float *)malloc(sizeof(float) * m);
    l_d = (float *)malloc(sizeof(float) * m);
    t = (float *)malloc(sizeof(float) * m * 2);
    tz = (float *)malloc(sizeof(float) * m);
    u_buff = (float *)malloc(sizeof(float) * lengthSeriesDataset);
    l_buff = (float *)malloc(sizeof(float) * lengthSeriesDataset);

    for (int i = 0; i < m; i++)
    {
        cb[i] = 0;
        cb1[i] = 0;
        cb2[i] = 0;
    }
    for (long int nS = startnumber; nS < stopnumber; nS++)
    {
        /// Read the buffer, which contains a single series in the dataset
        int k = 0;
        long long I, jj; /// the starting index of the data in current chunk of size EPOCH
        int kim = 0, keogh = 0, keogh2 = 0;
        //fread(buffer, lengthSeriesDataset*sizeof(float), 1, fp);
        //memcpy(buffer,&(rawfile[nS*lengthSeriesDataset]), sizeof(float)*lengthSeriesDataset);

        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop is end.

        lower_upper_lemire(&(rawfile[nS * lengthSeriesDataset]), lengthSeriesDataset, r, l_buff, u_buff);

        /// Do main task here..
        ex = 0;
        ex2 = 0;
        for (int i = 0; i < lengthSeriesDataset; i++)
        {
            /// A bunch of data has been read and pick one of them at a time to use
            d = (&(rawfile[nS * lengthSeriesDataset]))[i];

            /// Calcualte sum and sum square
            if (bNormalized)
            {
                ex += d;
                ex2 += d * d;
            }
            /// t is a circular array for keeping current data
            t[i % m] = d;

            /// Double the size for avoiding using modulo "%" operator
            t[(i % m) + m] = d;

            /// Start the task when there are more than m-1 points in the current chunk
            if (i >= m - 1)
            {
                //subsequences =subsequences+1;
                if (bNormalized)
                {
                    mean = ex / m;
                    std = ex2 / m;
                    std = sqrt(std - mean * mean);
                }
                /// compute the start location of the data in the current circular array, t
                jj = (i + 1) % m;
                /// the start location of the data in the current chunk
                I = i - (m - 1);

                /// Use a constant lower bound to prune the obvious subsequence
                lb_kim = lb_kim_hierarchy(t, q, jj, m, mean, std, bNormalized, bsf);

                if (lb_kim < bsf)
                {
                    /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                    /// uo, lo are envelop of the query.
                    lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, jj, m, mean, std, bNormalized, bsf);
                    if (lb_k < bsf)
                    {

                        if (bNormalized)
                        {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                            for (k = 0; k < m; k++)
                            {
                                tz[k] = (t[(k + jj)] - mean) / std;
                            }
                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelop for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I, u_buff + I, m, mean, std, bNormalized, bsf);
                        }
                        else
                        {
                            for (k = 0; k < m; k++)
                            {
                                tz[k] = t[(k + jj)];
                            }
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I, u_buff + I, m, mean, std, bNormalized, bsf);
                        }

                        if (lb_k2 < bsf)
                        {
                            /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                            /// Note that cb and cb2 will be cumulative summed here.
                            if (lb_k > lb_k2)
                            {
                                cb[m - 1] = cb1[m - 1];
                                for (k = m - 2; k >= 0; k--)
                                    cb[k] = cb[k + 1] + cb1[k];
                            }
                            else
                            {
                                cb[m - 1] = cb2[m - 1];
                                for (k = m - 2; k >= 0; k--)
                                    cb[k] = cb[k + 1] + cb2[k];
                            }

                            /// Compute DTW and early abandoning if possible

                            int pointChecked = 0;
                            dist = dtw(tz, q, cb, m, r, &pointChecked, bsf);
                            //totPointChecked = totPointChecked + pointChecked;

                            if (dist < bsf)
                            { /// Update bsf
                                /// loc is the real starting location of the nearest neighbor in the file
                                //omp_set_lock(&bsflock);
                                bsf = dist;
                                //omp_unset_lock(&bsflock);
                                //loc = (nS)*(lengthSeriesDataset-m+1) + i-m+1;
                            }
                        }
                    }
                }

                /// Reduce obsolute points from sum and sum square
                ex -= t[jj];
                ex2 -= t[jj] * t[jj];
            }
        }
    }
    printf("bsf is %f\n", bsf);
    free(cb);
    free(cb1);
    free(cb2);
    free(u_d);
    free(l_d);
    free(t);
    free(u_buff);
    free(tz);
    free(l_buff);

    *(((transdata *)transferdata)->bsf) = bsf;
}
/// Main Function
int main(int argc, char *argv[])
{
    FILE *fp;     /// data file pointer
    FILE *qp;     /// query file pointer
    float bsf;    /// best-so-far
    float *t, *q; /// data array and query array
    int *order;   ///new order of the query
    float *u, *l, *qo, *uo, *lo, *tz, *cb, *cb1, *cb2, *u_d, *l_d;

    float d;
    long long i, j, totPointsNotPruned;
    long long subsequences;
    long long totPointChecked;
    float ex, ex2, mean, std;
    int m = -1, r = -1;
    unsigned long lengthSeriesDataset = 0, numberSeriesDataset = 0;
    long long loc = 0;
    float t1, t2;
    float dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    float *buffer, *u_buff, *l_buff;
    Index *Q_tmp;
    _Bool bNormalized = true;
    float epsilon = -1;
    int threadnumber = 1, queriesnumber = 1;
    /// If not enough input, display an error.
    if (argc <= 3)
        error(4);

    /// read size of the query
    if (argc > 3)
        m = atol(argv[3]);

    /// read warping windows
    if (argc > 4)
    {
        r = atoi(argv[4]);
    }
    if (argc > 5)
    {
        numberSeriesDataset = atoi(argv[5]);
    }
    if (argc > 6)
    {
        lengthSeriesDataset = atoi(argv[6]);
    }

    if (argc > 7)
    {
        bNormalized = (atoi(argv[7]) == 1);
    }

    if (argc > 8)
    {
        queriesnumber = atof(argv[8]);
    }
    if (argc > 9)
    {
        threadnumber = atof(argv[9]);
    }
    if (argc > 10)
    {
        epsilon = atof(argv[10]);
    }
    qp = fopen(argv[2], "r");

    if (qp == NULL)
        error(2);

    /// malloc everything here
    q = (float *)malloc(sizeof(float) * m);
    if (q == NULL)
        error(1);
    qo = (float *)malloc(sizeof(float) * m);
    if (qo == NULL)
        error(1);
    uo = (float *)malloc(sizeof(float) * m);
    if (uo == NULL)
        error(1);
    lo = (float *)malloc(sizeof(float) * m);
    if (lo == NULL)
        error(1);

    order = (int *)malloc(sizeof(int) * m);
    if (order == NULL)
        error(1);

    Q_tmp = (Index *)malloc(sizeof(Index) * m);
    if (Q_tmp == NULL)
        error(1);

    u = (float *)malloc(sizeof(float) * m);
    if (u == NULL)
        error(1);

    l = (float *)malloc(sizeof(float) * m);
    if (l == NULL)
        error(1);

    cb = (float *)malloc(sizeof(float) * m);
    if (cb == NULL)
        error(1);

    cb1 = (float *)malloc(sizeof(float) * m);
    if (cb1 == NULL)
        error(1);

    cb2 = (float *)malloc(sizeof(float) * m);
    if (cb2 == NULL)
        error(1);

    u_d = (float *)malloc(sizeof(float) * m);
    if (u == NULL)
        error(1);

    l_d = (float *)malloc(sizeof(float) * m);
    if (l == NULL)
        error(1);

    t = (float *)malloc(sizeof(float) * m * 2);
    if (t == NULL)
        error(1);

    tz = (float *)malloc(sizeof(float) * m);
    if (tz == NULL)
        error(1);

    buffer = (float *)malloc(sizeof(float) * lengthSeriesDataset);
    if (buffer == NULL)
        error(1);

    u_buff = (float *)malloc(sizeof(float) * lengthSeriesDataset);
    if (u_buff == NULL)
        error(1);

    l_buff = (float *)malloc(sizeof(float) * lengthSeriesDataset);
    if (l_buff == NULL)
        error(1);

    float totTime = 0;
    float totKim = 0;
    float totKeogh = 0;
    float totKeogh2 = 0;
    float TotAbandoningPowerDTW = 0;

    int numberQuery = 0;
    float *rawfile = (float *)malloc(sizeof(float) * lengthSeriesDataset * numberSeriesDataset);
    fp = fopen(argv[1], "rb");
    if (fp == NULL)
    {
        error(2);
    }
    fseek(fp, 0L, SEEK_SET);
    unsigned long read_number = fread(rawfile, sizeof(float), lengthSeriesDataset * numberSeriesDataset, fp);
    char *line = NULL;
    size_t len = 0;
    int read;

    char nameFileQueryResult[255];
    char const *nnameQueryCompleteStat = "completeStatTRILLION_DTW";
    time_t rawtime;
    time(&rawtime);
    struct tm *readable;
    readable = localtime(&rawtime);
    sprintf(nameFileQueryResult, "QueryResult_TRILLION_DTW_%d", readable->tm_hour);
    FILE *queryAnswerFile;
    FILE *completeStatFile;

    for (int qsn = 0; qsn < queriesnumber; qsn++)
    {
        totPointChecked = 0;
        subsequences = 0;
        numberQuery++;
        /// start the clock of the single query
        t1 = clock();

        time_t tsss;
        tsss = time(NULL);
        /// Read query file
        if (epsilon >= 0)
        {
            bsf = epsilon;
        }
        else
        {
            bsf = FLT_MAX;
        }

        int totAnswerEpsilon = 0;
        long long posBestDist = 0;
        float bestDist = FLT_MAX;

        i = 0;
        j = 0;
        ex = ex2 = 0;
        fread(q, lengthSeriesDataset * sizeof(float), 1, qp);
        //char *olds = line;
        //char olddelim = ' ';
        //while(olddelim && *line) {
        //  while(*line && (' ' != *line)) line++;
        //  *line ^= olddelim = *line;
        //cb(olds);
        //   d = atof(olds);
        ////    ex += d;
        //  ex2 += d*d;
        //   q[i] = d;
        //   i++;
        //   *line++ ^= olddelim;
        //    olds = line;
        // }

        /// Do z-normalize the query, keep in same array, q
        if (bNormalized)
        {
            // mean = ex/m;
            //  std = ex2/m;
            //   std = sqrt(std-mean*mean);
            //   for( i = 0 ; i < m ; i++ )
            //       q[i] = (q[i] - mean)/std;
        }

        /// Create envelop of the query: lower envelop, l, and upper envelop, u
        lower_upper_lemire(q, m, r, l, u);

        /// Sort the query one time by abs(z-norm(q[i]))
        for (i = 0; i < m; i++)
        {
            Q_tmp[i].value = q[i];
            Q_tmp[i].index = i;
        }
        qsort(Q_tmp, m, sizeof(Index), comp);

        /// also create another arrays for keeping sorted envelop
        for (i = 0; i < m; i++)
        {
            int o = Q_tmp[i].index;
            order[i] = o;
            qo[i] = q[o];
            uo[i] = u[o];
            lo[i] = l[o];
        }

        /// Initial the cummulative lower bound
        for (i = 0; i < m; i++)
        {
            cb[i] = 0;
            cb1[i] = 0;
            cb2[i] = 0;
        }

        j = 0; /// the starting index of the data in the circular array, t
        ex = ex2 = 0;
        int k = 0;
        long long I; /// the starting index of the data in current chunk of size EPOCH
        int kim = 0, keogh = 0, keogh2 = 0;
        unsigned long nS;

        transdata transferdata[threadnumber];
        float bsfarray[threadnumber];
        for (int i = 0; i < threadnumber; i++)
        {
            transferdata[i].rawfile = rawfile;
            transferdata[i].q = q;
            transferdata[i].u = u;
            transferdata[i].l = l;
            transferdata[i].uo = uo;
            transferdata[i].lo = lo;
            transferdata[i].qo = qo;
            transferdata[i].cb = cb;
            transferdata[i].cb1 = cb1;
            transferdata[i].cb2 = cb2;
            transferdata[i].tz = tz;
            transferdata[i].t = t;
            transferdata[i].l_d = l_d;
            transferdata[i].u_d = u_d;
            transferdata[i].l_buff = l_buff;
            transferdata[i].u_buff = u_buff;
            transferdata[i].buffer = buffer;
            transferdata[i].bsf = &bsfarray[i];
            transferdata[i].m = m;
            transferdata[i].r = r;
            transferdata[i].order = order;
            transferdata[i].bNormalized = bNormalized;
            transferdata[i].startnumber = i * (numberSeriesDataset / threadnumber);
            transferdata[i].stopnumber = (i + 1) * (numberSeriesDataset / threadnumber);
        }

        //#pragma omp parallel for num_threads(3)
        pthread_t threadid[threadnumber];
        for (i = 0; i < threadnumber; i++)
        {
            pthread_create(&(threadid[i]), NULL, workerthread, (void *)&(transferdata[i]));
        }
        for (i = 0; i < threadnumber; i++)
        {
            pthread_join(threadid[i], NULL);
            if (bsfarray[i] < bsf)
                bsf = bsfarray[i];
        }
        bestDist = bsf;

        /*
        for(nS=0;nS<0;nS++)
        {
            /// Read the buffer, which contains a single series in the dataset

            //fread(buffer, lengthSeriesDataset*sizeof(float), 1, fp);
	        //printf("the ns is %ld \n",nS);
            memcpy(buffer,&(rawfile[nS*lengthSeriesDataset]), sizeof(float)*lengthSeriesDataset);

            /// Data are read in chunk of size EPOCH.
            /// When there is nothing to read, the loop is end.
         
            lower_upper_lemire(buffer, lengthSeriesDataset, r, l_buff, u_buff);

            /// Do main task here..
            ex=0; ex2=0;
            for(i=0; i<lengthSeriesDataset; i++)
            {
                /// A bunch of data has been read and pick one of them at a time to use
                d = buffer[i];

                /// Calcualte sum and sum square
                if(bNormalized)
                {
                    ex += d;
                    ex2 += d*d;
                }
                /// t is a circular array for keeping current data
                t[i%m] = d;

                /// Double the size for avoiding using modulo "%" operator
                t[(i%m)+m] = d;

                /// Start the task when there are more than m-1 points in the current chunk
                if( i >= m-1 )
                {
                    subsequences =subsequences+1;
                    if(bNormalized)
                    {
                        mean = ex/m;
                        std = ex2/m;
                        std = sqrt(std-mean*mean);
                    }
                    /// compute the start location of the data in the current circular array, t
                    j = (i+1)%m;
                    /// the start location of the data in the current chunk
                    I = i-(m-1);

                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, bNormalized, bsf);

                    if (lb_kim < bsf)
                    {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelop of the query.
                        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bNormalized, bsf);
                        if (lb_k < bsf)
                        {
                           
                            
                            if(bNormalized)
                            {
                                /// Take another linear time to compute z_normalization of t.
                                /// Note that for better optimization, this can merge to the previous function.
                                for(k=0;k<m;k++)
                                {   tz[k] = (t[(k+j)] - mean)/std;
                                }
                                /// Use another lb_keogh to prune
                                /// qo is the sorted query. tz is unsorted z_normalized data.
                                /// l_buff, u_buff are big envelop for all data in this chunk
                                lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bNormalized, bsf);
                            }
                            else
                            {
                                for(k=0;k<m;k++)
                                {   tz[k] = t[(k+j)];
                                }
                                lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bNormalized, bsf);
                            }
                
                            if (lb_k2 < bsf)
                            {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                if (lb_k > lb_k2)
                                {
                                    cb[m-1]=cb1[m-1];
                                    for(k=m-2; k>=0; k--)
                                        cb[k] = cb[k+1]+cb1[k];
                                }
                                else
                                {
                                    cb[m-1]=cb2[m-1];
                                    for(k=m-2; k>=0; k--)
                                        cb[k] = cb[k+1]+cb2[k];
                                }

                                /// Compute DTW and early abandoning if possible
                                
                                
                                int pointChecked = 0;
                                dist = dtw(tz, q, cb, m, r, &pointChecked, bsf);
                                totPointChecked = totPointChecked + pointChecked;
                                
                                if(epsilon > 0)
                                {
                                    if( dist < bsf )
                                    {   
                                        totAnswerEpsilon++;
                                        if(dist < bestDist)
                                        {
                                            posBestDist = (nS)*(lengthSeriesDataset-m+1) + i-m+1;
                                            bestDist = dist;
                                        }
                                    }
                                }
                                else
                                {
                                    if( dist < bsf )
                                    {   /// Update bsf
                                        /// loc is the real starting location of the nearest neighbor in the file
                                        bsf = dist;
                                        loc = (nS)*(lengthSeriesDataset-m+1) + i-m+1;
                                    }
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;

                    /// Reduce obsolute points from sum and sum square
                    ex -= t[j];
                    ex2 -= t[j]*t[j];
                }
            }
   
        }
        */
        totPointsNotPruned = (long long)(lengthSeriesDataset - m + 1) * m * (numberSeriesDataset) - ((keogh + keogh2 + kim) * m);
        i = (long long)(lengthSeriesDataset - m + 1) * (long long)m * (long long)numberSeriesDataset;

        t2 = clock();
        printf("\n");
        double timeQuery = time(NULL) - tsss;

        if (epsilon > 0)
        {
            /// Note that loc and i are long long.
            cout << "Location best answer : " << posBestDist << endl;
            cout << "Best Distance : " << bestDist << endl;
            cout << "Search space points : " << i << endl;
            cout << "Total answers: " << totAnswerEpsilon << endl;
            cout << "Total Execution Time : " << timeQuery << " sec" << endl;
        }
        else
        {
            /// Note that loc and i are long long.
            cout << "Location : " << loc << endl;
            cout << "Distance : " << (bsf) << endl;
            cout << "Search space points  : " << i << endl;
            cout << "Total Execution Time : " << timeQuery << " sec" << endl;
        }

        float prunedKim = ((float)kim / (float)subsequences) * 100;
        float prunedKeogh = ((float)keogh / (float)subsequences) * 100;
        float prunedKeogh2 = ((float)keogh2 / (float)subsequences) * 100;

        float abandoningPowerDTW = 100;

        if (totPointsNotPruned > 0)
        {
            float abandoningPowerDTW = 100 - (((float)totPointChecked / (float)totPointsNotPruned) * 100);
        }

        /// printf is just easier for formating ;)
        printf("\n");
        printf("Pruned by LB_Kim    : %6.2f%%\n", prunedKim);
        printf("Pruned by LB_Keogh  : %6.2f%%\n", prunedKeogh);
        printf("Pruned by LB_Keogh2 : %6.2f%%\n", prunedKeogh2);
        printf("DTW abandoning power : %6.2f%%\n", abandoningPowerDTW);

        totTime = totTime + timeQuery;
        totKim = totKim + prunedKim;
        totKeogh = totKeogh + prunedKeogh;
        totKeogh2 = totKeogh2 + prunedKeogh2;
        TotAbandoningPowerDTW = TotAbandoningPowerDTW + abandoningPowerDTW;
        // save result query
        //...
        queryAnswerFile = fopen(nameFileQueryResult, "a");
        fprintf(queryAnswerFile, "\n");

        if (epsilon > 0)
        {
            fprintf(queryAnswerFile, "Query number %d of length %d results: \n", numberQuery, m);
            fprintf(queryAnswerFile, "on datasets: %s\n", argv[1]);
            fprintf(queryAnswerFile, "Pos best answer: %lld\n", posBestDist);
            fprintf(queryAnswerFile, "Dist best answer: %f \n", (bestDist));
            fprintf(queryAnswerFile, "Total answers: %d \n", totAnswerEpsilon);
            fprintf(queryAnswerFile, "Pruned by LB_Kim    : %6.2f%%\n", prunedKim);
            fprintf(queryAnswerFile, "Pruned by LB_Keogh  : %6.2f%%\n", prunedKeogh);
            fprintf(queryAnswerFile, "Pruned by LB_Keogh2 : %6.2f%%\n", prunedKeogh2);
            fprintf(queryAnswerFile, "DTW abandoning power : %6.2f%%\n", abandoningPowerDTW);
            fprintf(queryAnswerFile, "Total time : %6.2f\n", timeQuery);
            fclose(queryAnswerFile);
        }
        else
        {
            fprintf(queryAnswerFile, "Query number %d of length %d results: \n", numberQuery, m);
            fprintf(queryAnswerFile, "on datasets: %s\n", argv[1]);
            fprintf(queryAnswerFile, "Pos 1-knn: %lld\n", loc);
            fprintf(queryAnswerFile, "Dist 1-knn: %f \n", (bsf));
            fprintf(queryAnswerFile, "Pruned by LB_Kim    : %6.2f%%\n", prunedKim);
            fprintf(queryAnswerFile, "Pruned by LB_Keogh  : %6.2f%%\n", prunedKeogh);
            fprintf(queryAnswerFile, "Pruned by LB_Keogh2 : %6.2f%%\n", prunedKeogh2);
            fprintf(queryAnswerFile, "DTW abandoning power : %6.2f%%\n", abandoningPowerDTW);
            fprintf(queryAnswerFile, "Total time    : %6.2f\n", timeQuery);
            fclose(queryAnswerFile);
        }
        line = NULL;
        len = 0;
    }

    // save stat
    //...
    float avgtotTime = (float)totTime / (float)numberQuery;
    float avgtotKim = (float)totKim / (float)numberQuery;
    float avgtotKeogh = (float)totKeogh / (float)numberQuery;
    float avgtotKeogh2 = (float)totKeogh2 / (float)numberQuery;
    float avgDTWAbPower = (float)TotAbandoningPowerDTW / (float)numberQuery;

    completeStatFile = fopen(nnameQueryCompleteStat, "a");
    fprintf(completeStatFile, "\n");

    if (epsilon > 0)
    {
        if (bNormalized)
        {
            fprintf(completeStatFile, "Average results of Z_NORMALIZED epslion (=%f)query length: %d (%s), warping window: %d \n", epsilon, m, argv[2], r);
            printf("Average results of Z_NORMALIZED epslion (=%f)query length: %d (%s), warping window: %d \n", epsilon, m, argv[2], r);
        }
        else
        {
            fprintf(completeStatFile, "Average results of NON Z_NORMALIZED epslion (=%f) query length: %d (%s), warping window: %d \n", epsilon, m, argv[2], r);
            printf("Average results of NON Z_NORMALIZED epslion (=%f) query length: %d (%s), warping window: %d \n", epsilon, m, argv[2], r);
        }
    }
    else
    {
        if (bNormalized)
        {
            fprintf(completeStatFile, "Average results of Z_NORMALIZED query length: %d (%s), warping window: %d \n", m, argv[2], r);
            printf("Average results of Z_NORMALIZED query length: %d (%s), warping window: %d \n", m, argv[2], r);
        }
        else
        {
            fprintf(completeStatFile, "Average results of NON Z_NORMALIZED query length: %d (%s), warping window: %d \n", m, argv[2], r);
            printf("Average results of NON Z_NORMALIZED query length: %d (%s), warping window: %d \n", m, argv[2], r);
        }
    }

    fprintf(completeStatFile, "on datasets: %s\n", argv[1]);
    fprintf(completeStatFile, "Total time    : %6.2f seconds \n", avgtotTime);
    fprintf(completeStatFile, "Pruned by LB_Kim  : %6.2f%%\n", avgtotKim);
    fprintf(completeStatFile, "Pruned by LB_Keogh : %6.2f%%\n", avgtotKeogh);
    fprintf(completeStatFile, "Pruned by LB_Keogh2 : %6.2f%%\n", avgtotKeogh2);
    fprintf(completeStatFile, "DTW Abandoning power     : %6.2f%%\n", avgDTWAbPower);

    printf("on datasets: %s\n", argv[1]);
    printf("Total time    : %6.2f seconds \n", avgtotTime);
    printf("Pruned by LB_Kim  : %6.2f%%\n", avgtotKim);
    printf("Pruned by LB_Keogh : %6.2f%%\n", avgtotKeogh);
    printf("Pruned by LB_Keogh2 : %6.2f%%\n", avgtotKeogh2);
    printf("DTW Abandoning power     : %6.2f%%\n", avgDTWAbPower);

    fclose(fp);
    fclose(completeStatFile);
    free(rawfile);
    free(q);
    free(u);
    free(l);
    free(uo);
    free(lo);
    free(qo);
    free(cb);
    free(cb1);
    free(cb2);
    free(tz);
    free(t);
    free(l_d);
    free(u_d);
    free(l_buff);
    free(u_buff);
    free(order);
    free(buffer);
    fclose(qp);
    free(Q_tmp);
    return 0;
}
