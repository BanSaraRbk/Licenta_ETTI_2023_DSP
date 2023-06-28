
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1D.h"
#include <TSystem.h>
#include <TROOT.h>
#include "TSpline.h"
// #include "header.h"
// std::vector<double> Trapezoid(std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
// {
// std::vector<double> result, p, s, y, d;

// result.push_back(signal_NB[0]);
// p.push_back(signal_NB[0]);
// s.push_back(signal_NB[0]);
// y.push_back(signal_NB[0]);

// for (int i = 1; i < signal_NB.size(); i++)
// {
// double term2 = signal_NB[i - TR_rise_time_K];
// double term3 = signal_NB[i - Trap_Delay_L];
// double term4 = signal_NB[i - Trap_Delay_L - TR_rise_time_K];

// double di = signal_NB[i] - term2 - term3 + term4 - signal_NB[i - 1] + term2 + term3 - term4;
// d.push_back(di);

// result.push_back(result.back() + di);
// p.push_back(p.back() + di);
// s.push_back(s.back() + p[i] + 500 * di);
// y.push_back(y.back() + result[i]);
// }

// return y;
// }

void Trapezoid(std::vector<double> *result, std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
{
	// std::vector<double> result;
	result->clear();

	result->push_back(signal_NB[0]);

	for (int i = 0; i < signal_NB.size(); i++)
	{
		double term2 = signal_NB[i - TR_rise_time_K];
		double term3 = signal_NB[i - Trap_Delay_L];
		double term4 = signal_NB[i - Trap_Delay_L - TR_rise_time_K];

		result->push_back(result->back() + signal_NB[i] - term2 - term3 + term4);
	}
}

std::vector<double> Baseline(const std::vector<double> signal)
{
	std::vector<double> signal_NB;
	signal_NB.clear();
	int baseline = std::accumulate(signal.begin(), signal.begin() + 200, 0.0) / 200;

	std::cout << "Baseline is " << baseline << std::endl;

	std::vector<double> x;
	for (int i = 0; i < signal.size(); i++)
	{
		signal_NB.push_back(signal[i] - baseline);
	}

	return signal_NB;
}

std::vector<double> Differentiate_Signal(const std::vector<double> signal)
{
	std::vector<double> signal_diff, yN;
	signal_diff.push_back(signal[0]);
	yN.push_back(signal[0]);
	std::vector<double> x;
	double kI = 0.02;
	double kD = 0.009;
	for (int i = 1; i < signal.size() - 1; i++)
	{

		yN.push_back((yN.back() + signal[i] - signal[i - 1]) / (1 + kD));
		signal_diff.push_back((signal_diff.back() + (1 + kI) * yN[i] - yN[i - 1]) / (1 + kD));
	}

	return signal_diff;
}
TPad *createThreeGraphsOnPad(const std::vector<double> &x, const std::vector<double> &first_signal, const std::vector<double> &second_signal, const std::vector<double> &third_signal)
{
	auto pad = new TPad("pad name", "pad title", 0, 0, 1, 1);
	pad->Divide(3, 1, 0.01, 0.01);
	pad->Draw();
	pad->cd(1);

	auto sig_graph_1 = new TGraph(x.size(), x.data(), first_signal.data());
	pad->cd(1);
	sig_graph_1->SetTitle(Form("Differential signal"));
	sig_graph_1->Draw("APL");

	auto sig_graph_2 = new TGraph(x.size(), x.data(), second_signal.data());
	sig_graph_2->SetMarkerColor(kBlue);
	sig_graph_2->SetMarkerStyle(kFullCircle);
	sig_graph_2->SetTitle(Form("Tapezoid CAEN"));
	pad->cd(2);
	sig_graph_2->Draw("AP*");

	auto sig_graph_3 = new TGraph(x.size(), x.data(), third_signal.data());
	sig_graph_3->SetMarkerColor(kRed);
	sig_graph_3->SetMarkerStyle(kFullSquare);
	sig_graph_3->SetTitle(Form("Tapezoid RAW"));
	pad->cd(3);
	sig_graph_3->Draw("AP*");

	return pad;
}

void PHA_git()
{

	TString filename = "test_000031.root";
	TString treename = "events";
	auto run_file = new TFile(filename, "READ");
	TTree *run_tree;
	auto loc_canv = new TCanvas("runs_from_DAQ.pdf", "c1");
	auto loc_pad = new TPad("pad name", "pad title", 0, 0, 1, 1);
	loc_pad->Divide(3, 1);
	loc_pad->Draw();

	run_tree = run_file->Get<TTree>(treename);

	Int_t channel;
	ULong_t timestamp;
	Short_t energy;

	vector<double> *waveform1 = {};
	vector<double> *waveform2 = {};
	UInt_t entries, i;

	Double_t TR_Flat_top_M = 500, TR_rise_time_K = 2500;
	double_t Trap_Delay_L = TR_Flat_top_M + TR_rise_time_K;

	run_tree->SetBranchAddress("channel", &channel);
	run_tree->SetBranchAddress("energy", &energy);
	run_tree->SetBranchAddress("timestamp", &timestamp);
	run_tree->SetBranchAddress("waveform1", &waveform1);
	run_tree->SetBranchAddress("waveform2", &waveform2);

	entries = run_tree->GetEntries();
	printf("entries = %d\n", entries);

	std::vector<double> x, raw_signal, trapezoid_ifin, analog_adder, signal_diff2, output_signal_adder_negative, output_signal_adder;

	for (int i = 0; i < 1; i++)
	{
		run_tree->GetEntry(i);
		cout << "ce canal e " << channel << endl;
		int bin = 0;
		for (int j = 0; j < waveform2->size() - 8; j++)
		{
			x.push_back(j);
			raw_signal.push_back(((*waveform2)[j] + (*waveform2)[j + 1] + (*waveform2)[j + 2] + (*waveform2)[j + 3]) / 4);
			trapezoid_ifin.push_back(((*waveform1)[j] + (*waveform1)[j + 1] + (*waveform1)[j + 2] + (*waveform1)[j + 3] + (*waveform1)[j + 4] + (*waveform1)[j + 5] + (*waveform1)[j + 6] + (*waveform1)[j + 7]) / 8);
		}
	}

	std::vector<double> raw_signal_NB = Baseline(raw_signal);

	std::vector<double> signal_diff = Differentiate_Signal(raw_signal_NB, 0.02, 0.009);

	for (int i = 0; i < raw_signal_NB.size(); ++i)
	{
		signal_diff2.push_back(raw_signal_NB[i] * (0.019));
	}

	for (int i = 0; i < signal_diff.size(); ++i)
	{
		analog_adder.push_back(signal_diff2[i] + signal_diff[i]);
	}

	Trapezoid(&output_signal_adder, analog_adder, TR_rise_time_K, Trap_Delay_L);

	for (int i = 0; i < signal_diff.size(); ++i)
	{
		output_signal_adder_negative.push_back(output_signal_adder[i] * (-1.0));
	}

	TPad *pad = createThreeGraphsOnPad(x, signal_diff, trapezoid_ifin, output_signal_adder_negative);
};
