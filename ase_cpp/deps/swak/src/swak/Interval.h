#ifndef SWAK_INTERVAL_H
#define SWAK_INTERVAL_H

/* 
 * An interval class with the end being exclusive.
 *
 * Interval and Interval32 are typedef'ed
 */
template <typename T>
class Interval_t
{
public:
  Interval_t()
  {
    my_start = my_end = 0;
  }

  Interval_t(T start_, T end_or_len, bool second_is_len = false)
  {
    my_start = start_;

    if (second_is_len)
      my_end = my_start + end_or_len;
    else
      my_end = end_or_len;
    
  }

  void SetLeft(T left_) { my_start = left_; }
  void SetRight(T right_) { my_end = right_; }
  void SetStart(T start_) {  my_start = start_; }
  void SetEnd(T end_) { my_end = end_; }
  void SetLen(T len_) { my_end = my_start + len_; }

  T left() const { return my_start; }
  T start() const { return my_start; }
  T end() const { return my_end; }
  T right() const { return my_end; }
  T len() const { return my_end - my_start; }
  T size() const { return my_end - my_start; }

  bool operator<(const Interval_t<T> &other) const { return my_end <= other.my_start; }
  bool operator<=(const Interval_t<T> &other) const { return my_end <= other.my_end; }
  bool operator>(const Interval_t<T> &other) const { return my_start >= other.my_end; }
  bool operator>=(const Interval_t<T> &other) const { return my_start >= other.my_start; }
  bool operator==(const Interval_t<T> &other) const { return (my_start == other.my_start) && (my_end == other.my_end); }
  bool operator!=(const Interval_t<T> &other) const { return (my_start != other.my_start) || (my_end != other.my_end); }
  
  Interval_t<T> GetIntersection(const Interval_t<T> &other) const
  {
    return Interval_t<T>(max(left(), other.left()), min(right(), other.right()), false);
  }

  bool Overlaps(const Interval_t<T> &other) const 
  {
    return left() < other.right() && other.left() < right();
  }

  Interval_t<T> GetUnion(const Interval_t<T> &other) const
  {
    return Interval_t<T>(min(left(), other.left()), max(right(), other.right()), false);
  }

  bool ContainsPoint(T p) const
  {
    return (my_start <= p) && (p < my_end);
  }

  bool Contains(const Interval_t<T> &other) const
  {
    return (other.left() >= left()) && (other.right() <= right());
  }

  bool StrictlyContains(const Interval_t<T> &other) const
  {
    return (other.left() > left()) && (other.right() < right());
  }

private:
  T my_start, my_end;
};

typedef Interval_t<int64> Interval;
typedef Interval_t<int> Interval32;

template <typename T>
ostream& operator<< (ostream &os, const Interval_t<T> &t)
{
  os << t.start() << "-" << t.end();
  return os;
}

#endif
