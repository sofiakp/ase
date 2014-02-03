#ifndef SWAK_DENSITY_ESTIMATOR
#define SWAK_DENSITY_ESTIMATOR

#include "Swak.h"

class DensityEstimator
{
public:
  vector<double> values;

  void Init(const vector<double> &unsorted_vals)
  {
    values = unsorted_vals;
    sort(values.begin(), values.end());
  }

  double AtLeast(double x)
  {
    vector<double>::iterator low_it = lower_bound(values.begin(), values.end(), x);

    double pos = int(low_it - values.begin());

    // Return what % of data points are >= this value
    return (values.size() - pos) / double(values.size()); 
  }
};

#endif
