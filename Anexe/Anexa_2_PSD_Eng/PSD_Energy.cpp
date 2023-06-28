
{

	TCanvas *canvas = new TCanvas("c", "Signals");
	canvas->Divide(2, 1, 0, 0);
	TFile *f_positron = TFile::Open("./SDataR_run_same_settings_delila_60Co.root");
	TH1D *h_output = new TH1D("h_output", "Energy(Channel0)", 10000, 0, 400000);
	TH1D *h_digi = new TH1D("h_output", "short", 10000, 0, 400000);

	TTree *t_data = (TTree *)f_positron->Get("Data_R");
	UShort_t t_channel, t_board, t_energy, t_energy_short;
	UInt_t t_flags;
	ULong64_t t_timestamp, t_timestamp2, t_timestamp_ch0, t_timestamp_ch1;
	TArrayS *samples;
	UInt_t entries, i, baseline, val;

	int pretrigger = 496 / 2;
	int gatepre = 100 / 2;

	UInt_t short_gate = 80 / 2, long_gate = 300 / 2, sample_signal = 128;
	UInt_t baseline_start = pretrigger - gatepre - sample_signal;

	Double_t sampleNumber_ch0[1000], rawSample[1000], sampleNumber_ch1[1000], rawSample_ch1[1000], timestamp[1000], sampleNumber[1000];

	entries = t_data->GetEntries();
	std::cout << "No of entries are" << entries << std ::endl;
	t_data->SetBranchAddress("Samples", &samples);
	t_data->SetBranchAddress("Timestamp", &t_timestamp);
	t_data->SetBranchAddress("Energy", &t_energy);
	t_data->SetBranchAddress("EnergyShort", &t_energy_short);
	t_data->SetBranchAddress("Channel", &t_channel);
	int iter = 0;
	int sum_long, sum_short;
	for (int iEntry = 0; iEntry < entries; ++iEntry)
	{
		int k = 0;

		t_data->GetEntry(iEntry);
		if (t_channel == 0)
		{
			for (int j = 0; j < samples->fN; j++)
			{
				sampleNumber[j] = j;
				rawSample[j] = samples->At(j);
			}
		}
		h_digi->Fill(t_energy);
		double baseline = 0;
		for (int k = 0; k < 100; k++)
		{
			baseline += rawSample[k];
		}
		baseline /= 100;

		sum_short = 0;

		for (int k = pretrigger - gatepre; k < short_gate + pretrigger - gatepre; k++)
		{

			sum_short += baseline - rawSample[k];
		}

		sum_long = 0;
		for (int k = pretrigger - gatepre; k < long_gate + pretrigger - gatepre; k++)
		{

			sum_long = sum_long + baseline - rawSample[k];
			if (sum_long > 7000)
			{
				h_output->Fill(sum_long);
			}
		}

		if (iter % 1000 == 0)
		{
			std::cout << "Ne aflam la evenimentul " << iter << endl;

			h_output->Draw();

			canvas->Update();
		}
		iter = iter + 1;
	}

	h_output->GetXaxis()->SetTitle("ADC Channels");
	h_output->GetYaxis()->SetTitle("Counts");
	h_output->Draw();
}
