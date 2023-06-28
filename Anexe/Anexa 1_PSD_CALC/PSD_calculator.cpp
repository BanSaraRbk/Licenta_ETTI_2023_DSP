
{

	TCanvas *canvas = new TCanvas("c", "Signals");
	canvas->Divide(2, 1, 0, 0);
	TFile *f_positron = TFile::Open("./SDataR_run_same_settings_delila_60Co.root");
	TH1D *h_output = new TH1D("h_output", "Pulse Shape Discriminator", 10000, -2, 2);

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
	double PSD_digitizor = 0, PSD_calculat = 0;
	TH2F *hcol21 = new TH2F("hcol21", "Option COLZ", 40, -4, 4, 40, -20, 20);
	for (int iEntry = 0; iEntry < 5000; ++iEntry)
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
		}

		PSD_digitizor = (t_energy * 1.0 - t_energy_short * 1.0) / t_energy * 1.0;
		PSD_calculat = (sum_long * 1.0 - sum_short * 1.0) / (sum_long * 1.0);

		h_output->Fill(PSD_calculat);
	}
	h_output->GetXaxis()->SetTitle("PSD");
	h_output->GetYaxis()->SetTitle("Counts");

	h_output->Draw();

	PSD_digitizor = (t_energy * 1.0 - t_energy_short * 1.0) / t_energy * 1.0;
	PSD_calculat = (sum_long * 1.0 - sum_short * 1.0) / (sum_long * 1.0);

	std::cout << "PSD calculat de catre digitizor este egal cu :" << PSD_digitizor << std::endl;
	std::cout << "PSD calculat de catre noi este egal cu :" << PSD_calculat << std::endl;

	std::cout << "Diferenta dintre cele doua este egala cu " << PSD_digitizor - PSD_calculat << std::endl;
	cout << "Tchannel e " << t_channel;
}
