/********************************************************************************
*
* Copyright (C) 2007-2010 Bauhaus-Universitaet Weimar
*
*********************************************************************************
*
*  module     : nurbscurve_impl.hpp
*
*  description:
*
********************************************************************************/
#ifndef GPUCAST_MATH_NURBSCURVE_IMPL_HPP
#define GPUCAST_MATH_NURBSCURVE_IMPL_HPP

namespace npr {
namespace gpucast { 
namespace math {

  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  nurbscurve<point_t>::nurbscurve()
    : _degree  (0),
      _points  ( ),
      _knots   ( )
  {}


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  nurbscurve<point_t>::nurbscurve(nurbscurve const& rhs)
  : _degree(rhs._degree),
    _points(rhs._points.begin(), rhs._points.end()),
    _knots(rhs._knots.begin(), rhs._knots.end())
  {}


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  nurbscurve<point_t>::~nurbscurve()
  {}


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::swap(nurbscurve& rhs)
  {
    std::swap(_degree, rhs._degree);
    std::swap(_points, rhs._points);
    std::swap(_knots, rhs._knots);
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::print(std::ostream& os) const
  {
    os << "nurbscurve: " << " \n" << "points : ";
    std::copy(_points.begin(), _points.end(), std::ostream_iterator<point_t>(os, ", "));
    os << "\n knots: ";
    std::copy(_knots.begin(), _knots.end(), std::ostream_iterator<typename point_t::value_type>(os, ", "));
    os << "\n";
  }


//////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::add(point_t const& point)
  {
    _points.push_back(point);
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::add_knot(typename nurbscurve::value_type knot)
  {
    _knots.push_back(knot);
    std::sort(_knots.begin(), _knots.end());
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::add_knot(typename nurbscurve::value_type knot, size_t m)
  {
    _knots.insert(_knots.end(), m, knot);
    std::sort(_knots.begin(), _knots.end());
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::clear()
  {
    _knots.clear();
    _points.clear();
    _degree = 0;
  }

  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline point_t nurbscurve<point_t>::evaluate(value_type const& t) const
  {
    auto N =  evaluate_basis_function(t);

    point_t enumerator;
    enumerator.weight(0);
    value_type denominator = 0;

    for (size_t i = 0; i != _points.size(); ++i) {
      auto tmp = _points[i];

      enumerator += N[i]  * tmp;
      denominator += N[i] ;
    }
     
    return enumerator / denominator;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline std::vector<typename nurbscurve<point_t>::value_type> nurbscurve<point_t>::evaluate_basis_function(typename nurbscurve<point_t>::value_type const& t) const
  {
    struct scalar_field_2d {
      scalar_field_2d(size_t rows, size_t cols) 
      {
        data.reserve(rows);
        for (size_t r = 0; r != rows; ++r) {
          data.emplace_back(std::vector<value_type>(cols, 0));
        }
      }

      std::vector<value_type> operator()(size_t r) const {
        std::vector<value_type> result;
        result.reserve(data.size());
        for (auto i : data) {
          result.emplace_back(i[r]);
        }
        return result;
      }

      value_type& operator()(size_t r, size_t c) {
        return data[r][c];
      }
      value_type const& operator()(size_t r, size_t c) const{
        return data[r][c];
      }
      std::vector<std::vector<value_type>> data;
    };

    scalar_field_2d N(_knots.size() - order(), order());

    for (size_t i = 0; i < _knots.size(); ++i) {
      if (i == get_knotspan(t)) {
        N(i, 0) = 1;
      }
    }

    for (size_t i = 1; i < order(); ++i) 
    {
      for (size_t j = 0; j < _knots.size() - order(); ++j)
      {
        auto t_i = _knots[j];
        auto t_ip = _knots[j + i];
        auto t_ip1 = _knots[j + i + 1];
        auto t_i1 = _knots[j + 1];

        if ((t_ip == t_i) && (t_ip1 == t_i1)) {
          N(j, i) = 0;
        } else {
          if (t_ip == t_i) {
            N(j, i) = ((t_ip1 - t) / (t_ip1 - t_i1)) * N(j + 1, i - 1);
          } else {
            if (t_ip1 == t_i1) {
              N(j, i) = ((t - t_i) / (t_ip - t_i)) * N(j, i - 1);
            }
            else {
              N(j, i) = ((t - t_i) / (t_ip - t_i)) * N(j, i - 1) + ((t_ip1 - t) / (t_ip1 - t_i1)) * N(j + 1, i - 1);
            }
          }
        }
      }


    }
    return N(order()-1);
  }

  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline size_t nurbscurve<point_t>::get_knotspan(typename nurbscurve<point_t>::value_type const& t) const
  {
    for (size_t i = 0; i < _knots.size(); ++i) {
      if (i + 1 < _knots.size()) {
        if (_knots[i] <= t && t <= _knots[i + 1]) {
          return i;
        }
      }
      else {
        // knotspan not found
        throw std::runtime_error("knotspan not found");
      }
    }
    throw std::runtime_error("knotspan not found");
  }

  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline size_t nurbscurve<point_t>::get_multiplicity(value_type const& t) const
  {
    auto m = 0;

    for (auto i : _knots.size()) {  //knots_ => _knots
      if (_knots[i] == t) {
        ++m;
      }
    }
    return m;
  }

  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline bool
  nurbscurve<point_t>::verify() const
  {
    return (_knots.size() == _points.size() + order());
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline size_t
  nurbscurve<point_t>::degree() const
  {
    return _degree;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline size_t
  nurbscurve<point_t>::order() const
  {
    return _degree + 1;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::degree(size_t deg)
  {
    _degree = deg;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline void
  nurbscurve<point_t>::normalize_knotvector()
  {
    value_type mn = _knots.front();
    value_type mx = _knots.back();

    std::vector<value_type> tmp;

    for (typename std::vector<value_type>::iterator i = _knots.begin(); i != _knots.end(); ++i) {
      tmp.emplace_back((*i - mn)/(mx - mn));
    }

    std::swap(tmp, _knots);
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  template<typename iterator_t>
  inline void
  nurbscurve<point_t>::set_knotvector(iterator_t begin, iterator_t end)
  {
    _knots.clear();
    _knots.reserve(end-begin);
    std::copy(begin, end, std::back_inserter(_knots));
    //normalize_knotvector();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  template<typename iterator_t>
  inline void
  nurbscurve<point_t>::set_points(iterator_t begin, iterator_t end)
  {
    _points.clear();
    _points.reserve(end-begin);
    std::copy(begin, end, std::back_inserter(_points));
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline std::vector<point_t> const&
  nurbscurve<point_t>::points() const
  {
    return _points;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline std::vector<typename point_t::value_type> const&
  nurbscurve<point_t>::knots() const
  {
    return _knots;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline typename nurbscurve<point_t>::point_iterator
  nurbscurve<point_t>::begin()
  {
    return _points.begin();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline typename nurbscurve<point_t>::point_iterator
  nurbscurve<point_t>::end()
  {
    return _points.end();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline typename nurbscurve<point_t>::const_point_iterator
    nurbscurve<point_t>::begin() const
  {
    return _points.begin();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline typename nurbscurve<point_t>::const_point_iterator
  nurbscurve<point_t>::end() const
  {
    return _points.end();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline point_t const&
  nurbscurve<point_t>::front() const
  {
    return _points.front();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline point_t const&
  nurbscurve<point_t>::back() const
  {
    return _points.back();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline point_t&
  nurbscurve<point_t>::operator[] ( std::size_t index )
  {
    assert(index < _points.size());
    return _points[index];
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline point_t const&
  nurbscurve<point_t>::operator[] ( std::size_t index ) const
  {
    assert(index < _points.size());
    return _points[index];
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline std::size_t
  nurbscurve<point_t>::size() const
  {
    return _points.size();
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline nurbscurve<point_t>&
  nurbscurve<point_t>::operator=(nurbscurve const& rhs)
  {
    nurbscurve tmp(rhs);
    swap(tmp);
    return *this;
  }


  //////////////////////////////////////////////////////////////////////////////
  template<typename point_t>
  inline std::ostream& operator<<(std::ostream& os, nurbscurve<point_t> const& nc) {
    nc.print(os);
    return os;
  }

} //namespace math 
} //namespace gpucast
} //namespace npr

#endif // GPUCAST_MATH_NURBSCURVE_IMPL_HPP
