// classes example
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
using namespace std;

class Missile {
  std::vector<double> x, y;
  std::vector<double> sorted;
  std::vector<int> best_indexs;
  public:
    void set_values (vector<double>,vector<double>);
    double find_med(vector<double>);
    double st_dev(vector<double>, double);
    void best_fit (vector<double>,vector<double>);
    void quick_sort (vector<double> &,int, int);
    void find_outliers (vector<double>,vector<double>,int);
    void printvectors() {for(int i=0;i<x.size();++i){ cout << x[i] << " " << y[i] << endl; }}
    void print_sorted() {for(int i=0;i<sorted.size();++i){ cout << sorted[i] << endl; }}
    
    struct distance_object{
      double distance;
      int index;
    };

};



void Missile::set_values (vector<double> a, vector<double>  b) {
  x = a;
  y = b;
}

void Missile::find_outliers(vector<double> a, vector<double>  b, int n){
  //Select rand 2 pts
  // compute line equation based on chosen pts,
  // calculate distance of remaining n-2 pts to best fitting line & record median error
  // retain current trail solution if median error is less than all previous trials
  double lowest_med = 1000;
  vector<int> best_inputs;
  vector<double> best_distances;
  vector<double> best_sorted;
  vector<distance_object> dist_table;
  for (int i=0; i<n;++i){
    double x1_tmp, y1_tmp, x2_tmp, y2_tmp, m, b_y0, sum_distance, total_avg,stdev, max_dist, grubbs = 0;
    vector<distance_object> tmp_dist_table;
    int rnd = (rand() % a.size());
     x1_tmp = a[rnd];
     y1_tmp = b[rnd];
    int rnd2 = (rand() % a.size());
    while (rnd2 == rnd){
      rnd2 = (rand() % a.size());
    }
    x2_tmp = a[rnd2];
    y2_tmp = b[rnd2];
    m = (y2_tmp - y1_tmp) / (x2_tmp - x1_tmp);
    b_y0 = y1_tmp - m * x1_tmp;
    
    // cout << "y=" << m << "x+" << b_y0 << " using x["<< rnd << "],y[" << rnd << "] and x[" << rnd2 << "],y[" << rnd2 <<  "]" << endl;
    for (int j=0;j<a.size();++j){
        // if (j == rnd || j == rnd2 )
          // continue; // this step was to remove calculations in case if the numbers selected are repeated, it ended up causing more trouble than it solved.
        double x0, y0, distance = 0;
        x0 = a[j];
        y0 = b[j];
        distance = ((abs( y0 - m*x0 - b_y0))/sqrt(pow(m,2)+1));
        sum_distance += distance;
        sorted.push_back(distance);
        // if (distance > max_dist){
        //   max_dist = distance;
        // }
        distance_object tmp_obj;
        tmp_obj.index = j;
        tmp_obj.distance = distance;
        tmp_dist_table.push_back(tmp_obj);
        // cout << "Distance to point " << x0 << " " << y0 << " is: " << distance << endl; 
    }
    quick_sort(sorted,0,sorted.size());
    double med = find_med(sorted);

    if (med < lowest_med){
      lowest_med = med;
      best_sorted = sorted;
      dist_table = tmp_dist_table;
      // print_sorted();
      // cout << "Median is: " << med << endl;
    }



    sorted.clear();
    sum_distance = 0;
    total_avg = 0;
    stdev = 0;
    grubbs = 0;
    
  }
  // cout << "best_sorted size: " << best_sorted.size() << endl;
  // cout << "lowest_med is: " << lowest_med << endl;
  // cout << "best inputs are: " << endl;
  // cout << "size of vec: " <<  a.size() << endl;
  int counter = 0;
  for (int z=0;z<best_sorted.size();++z){
      for(int y=0;y<dist_table.size();++y){
        if (best_sorted[z] == dist_table[y].distance && sorted[z] <= lowest_med ){
          counter++;
          // cout << a[dist_table[y].index] << " " << b[dist_table[y].index] << " with distance " << dist_table[y].distance <<  endl;
          // cout << "best are @ " << dist_table[y].index << endl;
          // cout << "dist here is: " << dist_table[y].distance << endl;
          best_indexs.push_back(dist_table[y].index);
        }
      }
  }
  // cout << "counter: "<<counter << endl;
}

void Missile::quick_sort(vector<double> & v, int left, int right) {
      int i = left, j = right;
      double tmp;
      double pivot = v[(left + right) / 2];
      /* partition */
      while (i <= j) {
            while (v[i] < pivot)
                  i++;
            while (v[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = v[i];
                  v[i] = v[j];
                  v[j] = tmp;
                  i++;
                  j--;
            }
      };
 
      /* recursion */
      if (left < j)
            quick_sort(v, left, j);
      if (i < right)
            quick_sort(v, i, right);
}

double Missile::find_med(vector<double> v)
{
    int n = v.size() / 2;
    //nth_element(v.begin(), v.begin()+n, v.end());
    // cout << "median indexs at index: " << n << endl;
    return v[n];
}


void Missile::best_fit(vector<double> q, vector<double>  w){
  double sum_x, sum_y, sum_xy_square, sum_xx, avg_x, avg_y, m, b_y0, y_res, sum_y_res, res, sum_res, r_sqr = 0;
  vector<double> a, b;
  for (int i=0;i<best_indexs.size();++i){
    a.push_back(q[best_indexs[i]]);
    b.push_back(w[best_indexs[i]]);
    // cout << "using " << best_indexs[i] << endl;
  }

  // cout << "going to use " << endl;
  // for (int i=0;i<a.size();++i){
  //   cout << a[i] << " " << b[i] << endl;
  // }

  // cout << "a.size(): " << a.size() << endl;

  for (int i=0;i<a.size();++i){
    sum_x += a[i];
    sum_y += b[i];
    sum_xy_square += a[i] * b[i];
    sum_xx += a[i] * a[i];
  }

  // cout << "Sum of x: " << sum_x << " Sum of y: " << sum_y << endl;
  // cout << "Sum of xy^2: " << sum_xy_square << " Sum of x^2: " << sum_xx << endl;

  avg_x = sum_x / a.size();
  avg_y = sum_y / b.size();

  // cout << "avg_x: " << avg_x << " avg_y: " << avg_y <<  endl;

  m = ((a.size()* sum_xy_square - sum_x * sum_y)/(a.size() * sum_xx - sum_x*sum_x));
  b_y0 = avg_y - m * avg_x;

  // cout << "Slope equals: " << m << " and y intercept equals: " << b_y0 << endl;
  cout << m << " " << b_y0 <<  endl;

  // for statistical fun
  for (int i=0; i< a.size();++i)
  {
    y_res = pow( b[i] - b_y0 - (m * a[i]) ,2);
    sum_y_res += y_res;
    res = pow(b[i] - avg_y,2);
    sum_res += res;
  }

  r_sqr = (sum_res - sum_y_res) / sum_res;

  // cout << "correlation coefficient(r) " << sqrt(r_sqr) << endl;

}

int main (int argc, char *argv[]) {
  srand ( time(NULL) );
    std::string fname = argv[1];
  std::vector<double> x;
  std::vector<double> y; 
  ifstream inp(fname.c_str(), ios::in | ios::binary);
  double tmp_x, tmp_y;
  while (inp >> tmp_x >> tmp_y){
    x.push_back(tmp_x);
    y.push_back(tmp_y);
  }

  Missile miss;
  miss.set_values(x,y);
  // 1-(1-(1-.30)^2)^5 should result in 0.9654974749 % 
  miss.find_outliers(x,y,5);
  miss.best_fit(x,y);

  return 0;
}
