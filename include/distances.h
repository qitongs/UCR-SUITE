//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_DISTANCES_H
#define UCR_SUITE_DISTANCES_H

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cumulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double * A, double * B, double * cb, int m, int r, double bsf = INF)
{

    double * cost;
    double * cost_prev;
    double * cost_tmp;
    int i, j, k;
    double x, y, z, min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double *) malloc(sizeof(double) * (2 * r + 1));
    for(k = 0; k < 2 * r + 1; k ++)    cost[k] = INF;

    cost_prev = (double *) malloc(sizeof(double) * (2 * r + 1));
    for(k = 0; k < 2 * r + 1; k ++)    cost_prev[k] = INF;

    for(i = 0; i < m; i ++)
    {
        k = max(0,r - i);
        min_cost = INF;

        for(j = max(0,i - r); j <= min(m - 1,i + r); j ++, k ++)
        {
            /// Initialize all row and column
            if((i == 0) && (j == 0))
            {
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i],B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i + r < m - 1 && min_cost + cb[i + r + 1] >= bsf)
        {   free(cost);
            free(cost_prev);
            return min_cost + cb[i + r + 1];
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

#endif //UCR_SUITE_DISTANCES_H
