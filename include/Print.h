#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <iterator>

struct Graph
{
	std::vector <double> x;	
	std::vector <double> Function;	
	std::vector <double> NoisedFunct; 
	std::vector <double> FiltFunct; 

	void print()
	{
		
		std::cout << "x,f(x) =";
		
		for (int i = 0; i < x.size(); ++i)
		{
			std::cout << "(" << x[i]<<";" << Function[i]<<") ";
		}

		std::cout << "\n\nx,f1(x) =";
		
		for (int i = 0; i < x.size(); ++i)
		{
			std::cout << "(" << x[i] << ";" << NoisedFunct[i] << ") ";
		}
		std::cout << "\n\nx,f2(x) =";
		
		for (int i = 0; i < x.size(); ++i)
		{
			std::cout << "(" << x[i] << ";" << FiltFunct[i] << ") ";
		}
	}
};
struct Mertics
{
	double omega = 0.0;
	double delta = 0.0;
	double J = 0.0;
};
struct Save
{
	double r = 0.0;
	double h = 0.0;
	double distance = 0.0;
	std::vector <double> alpha;
	Graph graphic;
	Mertics metrics;
	void print()
	{
		std::string temp;
		size_t reserve = (alpha.size() + (alpha.end() - alpha.begin() - 1u)) * 8u + 4u;
		temp = std::to_string(alpha.front());
		for (size_t index = 1u; index < alpha.size() + (alpha.end() - alpha.begin() - 1u); ++index)
		{
			if (index >= alpha.size())
			{
				temp = temp + ',' + std::to_string(alpha[alpha.size() - (index - (alpha.size() - 1u)) - 1u]);
				continue;
			}
			temp = temp + ',' + std::to_string(alpha[index]);
		}
		std::cout << std::left << std::setprecision(6) << "| " << std::setw(10) << h << " | " << std::setw(10) << distance << " |" << std::setw(reserve) << temp << "| " << std::setw(10) << metrics.omega << " | " << std::setw(10) << metrics.delta << " |\n";
	}
};
struct Result
{
	std::vector <Save> saves;
	Save best_save;

	void print()
	{
		size_t reserve = (saves.front().alpha.size() + (saves.front().alpha.end() - saves.front().alpha.begin() - 1u)) * 8u + 4u;
		std::cout << std::left << std::setprecision(6) << "| " << std::setw(10) << "h" << " | " << std::setw(10) << "dis" << " | " << std::setw(reserve) << "alpha" << " | " << std::setw(10) << "w" << " | " << std::setw(10) << "d" << " |\n";
		for (Save save : saves)
		{
			save.print();
		}
		std::cout << "\n";
		print_best();
		std::cout << "\n";
	}
	void print_best()
	{
		std::cout << std::left << std::setprecision(6) << "| " << std::setw(10) << "h*" << " | " << std::setw(10) << "J" << " | " << std::setw(10) << "w" << " | " << std::setw(10) << "d" << " |\n";
		std::cout << std::left << std::setprecision(6) << "| " << std::setw(10) << best_save.h << " | " << std::setw(10) << best_save.metrics.J << " | " << std::setw(10) << best_save.metrics.omega << " | " << std::setw(10) << best_save.metrics.delta << " |\n";
	}
};

struct Laba
{
	Result resultR3; 
	Result resultR5; 
};
