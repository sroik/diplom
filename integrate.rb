#!/usr/bin/ruby

# recursive factorial calculation
def exp(x=1)
  (Math.exp(x))
end

def factorial(n)
  (1..n).inject(:*) || 1.0
end

# exact exp_value
def exp_integral_value(t, m)
  e1 = exp(-2.0 * t)
  e2 = exp(-t*(m + 1.0))
  denominator = (m**2.0 - 1.0)
  (2.0*e2 - e1*(m + 1.0) + m - 1.0) / denominator
end

# accurate value
def exact_value(t, l, m)
  multiplier = -(l**2.0) / 4.0
  inside = multiplier * exp_integral_value(t, m)
  Math.exp(inside)
end

# under_sum value
def under_sum(t, lambda, m, j, k)
  numerator = lambda**(2.0*k)
  denominator = (2.0**k) * factorial(k) * ((m**2.0 - 1.0)**k)
  multiplier = numerator / denominator
  degree = 2**(j + 1.0)
  e1 = exp(-t / degree)
  e2 = exp(-m*t / degree)
  e3 = exp(-t*(m + 1.0) / degree)
  multiplier * ((e1 - e2)**k) * ((e3 - 1.0)**k)
end

# sum value
def sum_value(t, lambda, m, n, j)
  (0..n).inject(0) { |sum, k|
    sum + under_sum(t, lambda, m, j, k)
  }
end

def approximation_prod(t, lambda, m, n, l)
  (0...l).inject(1.0) do |prod, j|
    prod*sum_value(t, lambda, m, n, j)**(2.0**j)
  end
end

# approximation
# using prod formula
def approximation(t, lambda, m, n, l)
  prod = approximation_prod(t, lambda, m, n, l)
  exact_t = t / (2.0**l)
  exact = exact_value(exact_t, lambda, m)**(2.0**l)
  prod * exact
end

####
# tests running
####
def run_test(t, lambda, m, n, l)
  exact = exact_value(t, lambda, m)
  approx = approximation(t, lambda, m, n, l)
  puts "Test with T: #{t}, lambda: #{lambda}, m: #{m}, N: #{n}, l: #{l}"
  puts "exact value: #{exact}"
  puts "approximation: #{approx}"
end

t_values = [
  2,
  4,
  8,
  16
]

t_values.each do |t|
  run_test(t, 0.2, 2.0, 5, 5)
end

t_values.each do |t|
  run_test(t, 0.5, 10.0, 10, 10)
end

t_values.each do |t|
  run_test(t, 0.9, 15.0, 15, 15)
end
