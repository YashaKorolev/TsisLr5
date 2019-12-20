#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <vector>
#include <random>
#include <algorithm>
#include "Print.h"
#define Pi 3.141592653589793238463
struct Parameters
{
	std::pair<double, double> Span = { 0.0, Pi };
	size_t K = 100; 
	double Noise = 0.5; 
	double P = 0.95; 
	double eps = 0.01; 
	double J = 90000000.0;
	size_t L = 10u; 
	std::vector <size_t> Radius{ 3,5 }; 
	std::vector <size_t> M;
	std::vector <double> alpha;
	std::vector <double> xk;
	std::vector <double> Xk;
	std::vector <double> NoisedXk;
	std::vector <double> FiltXk;
	std::vector <double> Lambda; 
	std::pair<double, double> OmegaDelta{ 90000000.0,90000000.0 };

	void init()
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution <double> dis(-Noise / 2.0, Noise / 2.0);

		for (size_t index = 0u; index <= L; ++index)
		{
			Lambda.push_back(0.1 * index);
		}
		
		for (size_t index = 0u; index <= K; ++index)
		{
			xk.push_back(XkCalc(index));
		}
		
		for (size_t index = 0u; index <= K; ++index)
		{
			Xk.push_back(signal(xk[index]));
			NoisedXk.push_back(signal_plus_noise(xk[index], dis(gen)));
		}
		
		for (size_t index = 0u; index < Radius.size(); ++index)
		{
			M.push_back((Radius[index] - 1u) / 2u);
		}
		std::cout << "";
	}

	const double XkCalc(size_t k)
	{
		if (k >= 0u && k <= 100u)
		{
			return Span.first + k * (Span.second - Span.first) / K;
		}
		
	}
	
	const double signal(double Xk_)
	{
		return std::sin(Xk_) + 0.5;
	}
	
	const double signal_plus_noise(double Xk_, double noise)
	{
		return std::sin(Xk_) + 0.5 + noise;
	}
	
	const double getN()
	{
		return (log(1.0 - P) / log(1.0 - (eps / (Span.second - Span.first))));
	}
	
	double geom_sq(std::vector<double> alpha, size_t M, size_t k)
	{
		double sum = 0.0;
		if (k < M || k > K - M)
		{
			return Xk[k];
		}
		size_t degree = 0u;
		double AlphaTemp = 0.0;

		for (size_t index = k - M; index <= k + M; ++index)
		{
			degree = index + M - k;
			if (degree >= alpha.size())
			{
				degree = alpha[alpha.size() - (degree - (alpha.size() - 1u)) - 1u];
			}
			sum += pow(NoisedXk[index], 2) * alpha[degree];
		}
		return std::sqrt(sum);
	}
	double sum(std::vector<double> alpha, size_t a, size_t b)
	{
		double temp = 0.0;
		for (size_t index = a; index <= b; ++index)
		{
			if (index + 1u > alpha.size())
			{
				temp += alpha[index - alpha.size()];
				continue;
			}
			temp += alpha[index];
		}
		return temp;
	}
	void FiltFx(std::vector <double> input_alpha, size_t M)
	{
		FiltXk.clear();
		for (size_t index = 0u; index <= K; ++index)
		{
			FiltXk.push_back(geom_sq(input_alpha, M, index));
		}
	}
};

class job
{
	double N;
	Laba MyLaba;
	std::vector <Save> save;
	Parameters par;
public:
	job()
	{
		par.init();
		N = par.getN();
	}
	void pass()
	{
		Save temp_save;
		std::pair<double, double> TempOmegaDelta;
		double temp_J = 0.0;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution <double> dis(0.0, 1.0);
		std::vector<double> AlphaTemp;

		double temp;

		temp_save.graphic.x = par.xk;
		temp_save.graphic.Function = par.Xk;
		temp_save.graphic.NoisedFunct = par.NoisedXk;

		for (size_t jump = 0u; jump < par.Radius.size(); ++jump)
		{
			double r = par.Radius[jump];
			double M = par.M[jump];

			temp_save.r = par.Radius[jump];
			for (size_t l = 0u; l < par.Lambda.size(); ++l) 
			{
				temp_save.h = par.Lambda[l];
				for (size_t index = 0u; index < N; ++index)
				{
					AlphaTemp.reserve(M + 1u);
					AlphaTemp.resize(M + 1u);
					AlphaTemp.back() = dis(gen);
					if (par.M[jump] >= 2u)
					{
						for (size_t count = 1u; count < par.M[jump]; ++count)
						{
							std::uniform_real_distribution <double> temp_dis(0.0, 1.0 - par.sum(AlphaTemp, M, r - M - 1));
							AlphaTemp.at(count) = 0.5 * temp_dis(gen);
						}
					}
					AlphaTemp.front() = 0.5 * (1.0 - par.sum(AlphaTemp, 1u, r - 2u));
					
					par.FiltFx(AlphaTemp, M);
					
					double summ = 0.0;
					for (size_t kol = 1u; kol <= par.K; ++kol)
					{
						summ += abs(par.FiltXk[kol] - par.FiltXk[kol - 1u]);
					}
					TempOmegaDelta.first = summ; summ = 0.0;
					
					for (size_t kol = 0u; kol <= par.K; ++kol)
					{
						summ += abs((par.FiltXk[kol] - par.NoisedXk[kol]));
					}
					TempOmegaDelta.second = summ / par.K;

					temp_J = (par.Lambda[l] * TempOmegaDelta.first) + (1 - par.Lambda[l]) * TempOmegaDelta.second;

					if (temp_J < par.J)
					{
						par.J = temp_J;
						par.OmegaDelta = TempOmegaDelta;
						temp_save.alpha = AlphaTemp;
						temp_save.metrics.delta = TempOmegaDelta.second;
						temp_save.metrics.omega = TempOmegaDelta.first;
						temp_save.metrics.J = temp_J;
						temp_save.graphic.FiltFunct = par.FiltXk;
						temp_save.distance = abs(temp_save.metrics.delta) + abs(temp_save.metrics.omega);
					}
				}
				par.J = 90000000.0;
				par.OmegaDelta = { 90000000.0,90000000.0 };
				save.push_back(temp_save);
			}

			if (jump == 0u)
			{
				MyLaba.resultR3.saves = save;
				save.clear();
			}
			else
			{
				MyLaba.resultR5.saves = save;
				save.clear();
			}
		}

		find_best();
		std::cout << "\t\tR = 3 \n";
		MyLaba.resultR3.print();
		std::cout << "\t\tR = 5 \n";
		MyLaba.resultR5.print();
		std::cout << "\n\n For R = 3\n";
		MyLaba.resultR3.best_save.graphic.print();
		std::cout << "\n\n\nFor R = 5 \n";
		MyLaba.resultR5.best_save.graphic.print();
	}
	void find_best()
	{
		Save temp_best;
		double best_distance = 999999900.0;
		for (Save save : MyLaba.resultR3.saves)
		{
			if (save.distance < best_distance)
			{
				best_distance = save.distance;
				temp_best = save;
			}
		}
		MyLaba.resultR3.best_save = temp_best;
		best_distance = 999999900.0;
		for (Save save : MyLaba.resultR5.saves)
		{
			if (save.distance < best_distance)
			{
				best_distance = save.distance;
				temp_best = save;
			}
		}
		MyLaba.resultR5.best_save = temp_best;
	}
};
