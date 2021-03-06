/********************************************************************
 ********************************************************************
 ** C++ class implementation of the Hungarian algorithm by David Schwarz, 2012
 **
 **
 ** O(n^3) implementation derived from libhungarian by Cyrill Stachniss, 2004
 **
 **
 ** Solving the Minimum Assignment Problem using the 
 ** Hungarian Method.
 **
 ** ** This file may be freely copied and distributed! **
 **
 **
 ** This file is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied 
 ** warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 ** PURPOSE.  
 **
 ********************************************************************
 ********************************************************************/

#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>

#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#define verbose (1)


using std::vector;

typedef enum {
	HUNGARIAN_MODE_MINIMIZE_COST,
	HUNGARIAN_MODE_MAXIMIZE_UTIL,
} MODE;

typedef enum {
	HUNGARIAN_NOT_ASSIGNED,
	HUNGARIAN_ASSIGNED,
} ASSIGN;


template<class T>
class Hungarian
{

public:
	/** This method initialize the hungarian_problem structure and init 
	 *  the  cost matrices (missing lines or columns are filled with 0).
	 *  It returns the size of the quadratic(!) assignment matrix. **/

	Hungarian();
	Hungarian(const vector<vector<T> >&, int, int, MODE);

	int init(const vector<vector<T> >& input_matrix,
			   int rows, 
			   int cols, 
			   MODE mode);

	/** This method computes the optimal assignment. **/
	bool solve();

	/** Accessor for the cost **/
	int cost() const;

	/** Reference accessor for assignment **/
	const vector<vector<T> >& assignment() const;

	/** Print the computed optimal assignment. **/
	void print_assignment();

	/** Print the cost matrix. **/
	void print_cost();

	/** Print cost matrix and assignment matrix. **/
	void print_status();

protected:
	bool check_solution(const vector<T>& row_dec, const vector<T>& col_inc, const vector<int>& col_vertex);
	bool assign_solution(const vector<T>& row_dec, const vector<T>& col_inc, const vector<int>& col_vertex);

private:
	T m_epsilon;

	T m_cost;
	int m_rows;
	int m_cols;
	vector<vector<T> > m_costmatrix;
	vector<vector<int> > m_assignment;

};




template<class T>
Hungarian<T>::Hungarian()
{
	m_epsilon = 1e-9;

	//much ado about nothing
	m_rows = 1;
	m_cols = 1;

	m_costmatrix.resize(m_rows, vector<T>(m_cols,0));
	m_assignment.resize(m_rows, vector<int>(m_cols,0));
}

template<class T>
Hungarian<T>::Hungarian(const vector<vector<T> >& input_matrix, int rows, int cols, MODE mode)
{
	m_epsilon = 1e-9;

	int i,j, org_cols, org_rows;
	T max_cost;
	max_cost = 0;

	org_cols = cols;
	org_rows = rows;

	// is the matrix square?
	// if no, expand with 0-cols / 0-cols

	if(rows!=cols)
	{
		rows = std::max(cols, rows);
		cols = rows;
	}

	m_rows = rows;
	m_cols = cols;

	m_costmatrix.resize(rows, vector<T>(cols,0));
	m_assignment.resize(rows, vector<int>(cols,0));

	for(i=0; i<m_rows; i++)
	{
		for(j=0; j<m_cols; j++)
		{
			m_costmatrix[i][j] =  (i < org_rows && j < org_cols) ? input_matrix[i][j] : 0;
			m_assignment[i][j] = 0;

			if (max_cost < m_costmatrix[i][j])
			{
				max_cost = m_costmatrix[i][j];
			}
		}
	}


	if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL)
	{
		for(i=0; i<m_rows; i++)
		{
			for(j=0; j<m_cols; j++)
			{
				m_costmatrix[i][j] =  max_cost - m_costmatrix[i][j];
			}
		}
	}
	else if (mode == HUNGARIAN_MODE_MINIMIZE_COST)
	{
		// nothing to do
	}
	else
		fprintf(stderr,"%s: unknown mode. Mode was set to HUNGARIAN_MODE_MINIMIZE_COST !\n", __FUNCTION__);
}


template<class T>
void hungarian_print_matrix(const vector<vector<T> >& C, int rows, int cols)
{
	int i,j;
	fprintf(stderr , "\n");
	for(i=0; i<rows; i++)
	{
		fprintf(stderr, " [");
		for(j=0; j<cols; j++)
		{
			fprintf(stderr, "%5f ",C[i][j]);
		}
		fprintf(stderr, "]\n");
	}
	fprintf(stderr, "\n");
}

template<class T>
void Hungarian<T>::print_assignment() {
	hungarian_print_matrix(m_assignment, m_rows, m_cols) ;
}

template<class T>
void Hungarian<T>::print_cost() {
	hungarian_print_matrix(m_costmatrix, m_rows, m_cols) ;
}

template<class T>
void Hungarian<T>::print_status()
{

	fprintf(stderr,"cost:\n");
	print_cost();

	fprintf(stderr,"assignment:\n");
	print_assignment();

}

template<class T>
int Hungarian<T>::init(const vector<vector<T> >& input_matrix, int rows, int cols, MODE mode)
{

	int i,j, org_cols, org_rows;
	T max_cost;
	max_cost = 0;

	org_cols = cols;
	org_rows = rows;

	// is the number of cols  not equal to number of rows ?
	// if yes, expand with 0-cols / 0-cols
	rows = std::max(cols, rows);
	cols = rows;

	m_rows = rows;
	m_cols = cols;

	m_costmatrix.resize(rows, vector<T>(cols,0));
	m_assignment.resize(rows, vector<int>(cols,0));

	for(i=0; i<m_rows; i++)
	{
		for(j=0; j<m_cols; j++)
		{
			m_costmatrix[i][j] =  (i < org_rows && j < org_cols) ? input_matrix[i][j] : 0;
			m_assignment[i][j] = 0;

			if (max_cost < m_costmatrix[i][j])
				max_cost = m_costmatrix[i][j];
		}
	}


	if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL) {
		for(i=0; i<m_rows; i++) {
			for(j=0; j<m_cols; j++) {
				m_costmatrix[i][j] =  max_cost - m_costmatrix[i][j];
			}
		}
	}
	else if (mode == HUNGARIAN_MODE_MINIMIZE_COST) {
		// nothing to do
	}
	else
		fprintf(stderr,"%s: unknown mode. Mode was set to HUNGARIAN_MODE_MINIMIZE_COST !\n", __FUNCTION__);

	return rows;
}

template<class T>
bool Hungarian<T>::check_solution(const vector<T>& row_dec, const vector<T>& col_inc, const vector<int>& col_vertex)
{
	int k, l, m, n;

	m = m_rows;
	n = m_cols;
	// Begin doublecheck the solution 23
	for (k=0;k<m;k++)
		for (l=0;l<n;l++)
			if (m_costmatrix[k][l]<row_dec[k]-col_inc[l])
				return false;

	for (k=0;k<m;k++)
	{
		l=col_vertex[k];
		if (l<0 || fabs(m_costmatrix[k][l]-(row_dec[k]-col_inc[l]))>m_epsilon)
			return false;
	}
	k=0;

	for (l=0;l<n;l++)
	{
		if (col_inc[l])
		{
			k++;
		}
	}
	if (k>m)
	{
		return false;
	}

	//everything checks out, then
	return true;
	// End doublecheck the solution 23
}

template<class T>
bool Hungarian<T>::assign_solution(const vector<T>& row_dec,const vector<T>&  col_inc, const vector<int>&  col_vertex)
{
	// End Hungarian algorithm 18
	int i, j, k, l, m, n;

	m = m_rows;
	n = m_cols;

	for (i=0;i<m;++i)
	{
		m_assignment[i][col_vertex[i]]=HUNGARIAN_ASSIGNED;
		/*TRACE("%d - %d\n", i, col_vertex[i]);*/
	}
	for (k=0;k<m;++k)
	{
		for (l=0;l<n;++l)
		{
			/*TRACE("%d ",m_costmatrix[k][l]-row_dec[k]+col_inc[l]);*/
			m_costmatrix[k][l]=m_costmatrix[k][l]-row_dec[k]+col_inc[l];
		}
		/*TRACE("\n");*/
	}
	for (i=0;i<m;i++)
	{
		m_cost+=row_dec[i];
	}
	for (i=0;i<n;i++)
	{
		m_cost-=col_inc[i];
	}
	if (verbose)
		fprintf(stderr, "Cost is %f\n",m_cost);

	return true;

}

template<class T>
bool Hungarian<T>::solve()
{
	int i, j, m, n, k, l, t, q, unmatched;
	T s, cost;

	m = m_rows;
	n = m_cols;

	T INF = std::numeric_limits<T>::max();

	//vertex alternating paths,
	vector<int> col_vertex(m), row_vertex(n), unchosen_row(m), parent_row(n),
					slack_row(m);
	vector<T>   slack(n),  col_inc(n), row_dec(m);

	cost=0;

	for (i=0;i<m_rows;i++)
	{
		col_vertex[i]=0;
		unchosen_row[i]=0;
		row_dec[i]=0;
		slack_row[i]=0;
	}

	for (j=0;j<m_cols;j++)
	{
		row_vertex[j]=0;
		parent_row[j] = 0;
		col_inc[j]=0;
		slack[j]=0;
	}

	//Double check assignment matrix is 0
	m_assignment.assign(m, vector<int>(n, HUNGARIAN_NOT_ASSIGNED));

	// Begin subtract column minima in order to start with lots of zeroes 12
	if (verbose)
	{
		fprintf(stderr, "Using heuristic\n");
	}

	for (l=0;l<n;l++)
	{
		s = m_costmatrix[0][l];

		for (k=1;k<m;k++)
		{
			if (m_costmatrix[k][l] < s)
			{
				s=m_costmatrix[k][l];
			}
			cost += s;
		}

		if (fabs(s)>=m_epsilon)
		{
			for (k=0;k<m;k++)
			{
				m_costmatrix[k][l]-=s;
			}
		}

		//pre-initialize state 16
		row_vertex[l]= -1;
		parent_row[l]= -1;
		col_inc[l]=0;
		slack[l]=INF;
	}
	// End subtract column minima in order to start with lots of zeroes 12

	// Begin initial state 16
	t=0;

	for (k=0;k<m;k++)
	{
		bool row_done = false;
		s=m_costmatrix[k][0];

		for (l=0;l<n;l++)
		{

			if(l > 0)
			{
				if (m_costmatrix[k][l] < s)
				{
					s = m_costmatrix[k][l];
				}
				row_dec[k]=s;
			}

//			if (fabs(s-m_costmatrix[k][l])>m_epsilon && row_vertex[l]<0)
			if (s == m_costmatrix[k][l] && row_vertex[l]<0)
			{
				col_vertex[k]=l;
				row_vertex[l]=k;

				if (verbose)
				{
					fprintf(stderr, "matching col %d==row %d\n",l,k);
				}
				row_done = true;
				break;
			}
		}

		if(!row_done)
		{
			col_vertex[k]= -1;

			if (verbose)
			{
				fprintf(stderr, "node %d: unmatched row %d\n",t,k);
			}

			unchosen_row[t++]=k;
		}

	}
	// End initial state 16

	bool checked = false;

	// Begin Hungarian algorithm 18

	//is matching already complete?
	if (t == 0)
	{
		checked = check_solution(row_dec, col_inc, col_vertex);
		if (checked)
		{
			//finish assignment, wrap up and done.
			bool assign = assign_solution(row_dec, col_inc, col_vertex);
			return true;
		}
		else
		{
			if(verbose)
			{
				fprintf(stderr, "Could not solve. Error.\n");
			}
			return false;
		}
	}

	unmatched=t;


	while (1)
	{
		if (verbose)
		{
			fprintf(stderr, "Matched %d rows.\n",m-t);
		}
		q=0;
		bool try_matching;
		while (1)
		{
			while (q<t)
			{
				// Begin explore node q of the forest 19

				k=unchosen_row[q];
				s=row_dec[k];
				for (l=0;l<n;l++)
				{
					if (fabs(slack[l])>=m_epsilon)
					{
						T del;
						del=m_costmatrix[k][l]-s+col_inc[l];
						if (del<slack[l])
						{
							if (fabs(del)<m_epsilon)
							{
								if (row_vertex[l]<0)
								{
									goto breakthru;
								}
								slack[l]=0;
								parent_row[l]=k;
								if (verbose)
									fprintf(stderr, "node %d: row %d==col %d--row %d\n",
											t,row_vertex[l],l,k);
								unchosen_row[t++]=row_vertex[l];
							}
							else
							{
								slack[l]=del;
								slack_row[l]=k;
							}
						}
					}
				}
				// End explore node q of the forest 19
				q++;
			}

			// Begin introduce a new zero into the matrix 21
			s=INF;
			for (l=0;l<n;l++)
			{
				if (fabs(slack[l])>m_epsilon && fabs(slack[l]-s)>m_epsilon)
				{
					s=slack[l];
				}
			}
			for (q=0;q<t;q++)
			{
				row_dec[unchosen_row[q]]+=s;
			}
			for (l=0;l<n;l++)
			{
				//check slack
				if (fabs(slack[l])>m_epsilon)
				{
					slack[l]-=s;
					if (fabs(slack[l])<=m_epsilon)
					{
						// Begin look at a new zero 22
						k=slack_row[l];
						if (verbose)
						{
							fprintf(stderr,
									"Decreasing uncovered elements by %f produces zero at [%d,%d]\n",
									s,k,l);
						}
						if (row_vertex[l]<0)
						{
							for (j=l+1;j<n;j++)
								if (fabs(slack[j])<=m_epsilon)
								{
									col_inc[j]+=s;
								}

							goto breakthru;
						}
						else
						{
							parent_row[l]=k;
							if (verbose)
							{
								fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_vertex[l],l,k);
								unchosen_row[t++]=row_vertex[l];
							}
						}
						// End look at a new zero 22
					}
				}
				else
				{
					col_inc[l]+=s;
				}
			}
			// End introduce a new zero into the matrix 21
		}

	  breakthru:
		// Begin update the matching 20
		if (verbose)
		{
			fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
		}
		while (1)
		{
			j=col_vertex[k];
			col_vertex[k]=l;
			row_vertex[l]=k;
			if (verbose)
			{
				fprintf(stderr, "rematching col %d==row %d\n",l,k);
			}
			if (j<0)
			{
				break;
			}
			k=parent_row[j];
			l=j;
		}
		// End update the matching 20
		if (--unmatched == 0)
		{
			checked = check_solution(row_dec, col_inc, col_vertex);
			if (checked)
			{
				//finish assignment, wrap up and done.
				bool assign = assign_solution(row_dec, col_inc, col_vertex);
				return true;
			}
			else
			{
				if(verbose)
				{
					fprintf(stderr, "Could not solve. Error.\n");
				}
				return false;
			}
		}

		// Begin get ready for another stage 17
		t=0;
		for (l=0;l<n;l++)
		{
			parent_row[l]= -1;
			slack[l]=INF;
		}
		for (k=0;k<m;k++)
		{
			if (col_vertex[k]<0)
			{
				if (verbose)
					fprintf(stderr, "node %d: unmatched row %d\n",t,k);
				unchosen_row[t++]=k;
			}
		}
		// End get ready for another stage 17
	}// back to while loop


}

//ACCESSORS

template<class T>
int Hungarian<T>::cost() const
{
	return m_cost;
}

template<class T>
const vector<vector<T> >& Hungarian<T>::assignment() const
{
	return m_assignment;
}


#endif



