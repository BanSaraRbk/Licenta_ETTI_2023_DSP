{

	TFile *f_positron = TFile::Open("./SDataR_run_same_settings_delila_60Co.root");

	TTree *t_data = (TTree *)f_positron->Get("Data_R");
	UShort_t t_channel, t_channel2, t_board, t_energy, t_energy2;
	UInt_t t_flags;
	ULong64_t t_timestamp, t_timestamp2, t_timestamp_ch0, t_timestamp_ch1;
	TArrayS *samples;
	UInt_t entries, i, baseline, val, j, k;
	;
	Double_t FineTS, fraction_cfd = 0.75;
	UInt_t delay_cfd = 3;
	int pretrigger = 1000 / 4;
	int gatepre = 48 / 4;

	UInt_t short_gate = 80 / 4, long_gate = 300 / 4, sample_signal = 128;
	UInt_t baseline_start = pretrigger - gatepre - sample_signal;

	TH1D *h_digi_output = new TH1D("h_digi_output", "Time Difference CAEN", 1000, -12000, 12000);
	TH1D *h_myoutput = new TH1D("h_myoutput", "short", 1000, 6000, 12000);

	h_digi_output->SetFillColor(kBlue);
	h_myoutput->SetFillColor(kBlue);

	TH1D *h_digi_FTS = new TH1D("h_digi_FTS", "Timestamp CAEN ", 200, -1000, 3000);
	TH1D *h_my_FTS = new TH1D("h_my_FTS", "h_my_FTS", 200, -1000, 3000);
	h_digi_FTS->SetFillColor(kBlue);
	h_my_FTS->SetFillColor(kBlue);

	auto loc_canv = new TCanvas("runs_from_delila", "c1", 1000, 800);
	loc_canv->Divide(2);
	Double_t sampleNumber_ch0[1000], rawSample_ch0[1000], sampleNumber_ch1[1000], rawSample_ch1[1000], timestamp[1000];
	std::vector<double> sampleNumber;
	entries = t_data->GetEntries();
	std::cout << "No of entries are" << entries << std ::endl;
	t_data->SetBranchAddress("Samples", &samples);
	t_data->SetBranchAddress("Timestamp", &t_timestamp);
	t_data->SetBranchAddress("Energy", &t_energy);
	t_data->SetBranchAddress("Channel", &t_channel);

	Double_t digitmp_0, digitmp_1;

	TGraph *ch0, *pulse0;

	int iter = 0;
	for (int iEntry = 0; iEntry < entries; iEntry++)
	{

		t_data->GetEntry(iEntry);

		if (t_channel == 0)
		{
			digitmp_0 = t_timestamp;
			h_digi_FTS->Fill((t_timestamp % 2000));
		}

		if (t_channel == 1)
		{
			digitmp_1 = t_timestamp;
		}

		h_digi_output->Fill(digitmp_1 - digitmp_0);

		if (iter % 1000 == 0)
		{
			std::cout << "Ne aflam la evenimentul " << iter << endl;

			loc_canv->cd(1);
			h_digi_output->Draw();
			loc_canv->cd(2);
			h_digi_FTS->Draw();

			loc_canv->Update();
		}
		iter = iter + 1;
	}
}