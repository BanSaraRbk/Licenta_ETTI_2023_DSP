{
	double Interpolation_orig(const Double_t *cfd, std::vector<double> sampleNumber, Double_t *fts)
	{
		bool flag = false;
		int a;
		double x0, y0, x1, y1, dx, dy, intercept, slope, ZC, Tfine;

		for (a = 0; a < sampleNumber.size() - 1; a++)
		{
			if (cfd[a] < -100)
			{
				flag = true;
			}

			if (flag == true && (cfd[a + 1] > 0))
			{
				break;
			}
		}

		x0 = a;
		y0 = cfd[a];
		x1 = a + 1;
		y1 = cfd[a + 1];
		dx = x1 - x0;
		dy = y1 - y0;
		slope = dy / dx;
		intercept = y1 - slope * x1;

		ZC = -intercept / slope;
		Tfine = ZC - x0;
		*fts = Tfine;
		return ZC;
	}

	TFile *f_positron = TFile::Open("./run128_0_photofission.root ");

	TTree *t_data = (TTree *)f_positron->Get("ELIADE_Tree");
	UChar_t t_channel, t_channel2, t_board, t_energy, t_energy2;
	UInt_t t_flags;
	ULong64_t t_timestamp, t_timestamp2, t_timestamp_ch0, t_timestamp_ch1;
	UShort_t data[10000];
	UInt_t recLength;
	UInt_t entries, i, baseline, val, j, k;
	;
	Double_t FineTS, fraction_cfd = 0.75;
	UInt_t delay_cfd = 12;
	int pretrigger = 1000 / 4;
	int gatepre = 48 / 4;

	UInt_t short_gate = 80 / 4, long_gate = 300 / 4, sample_signal = 128;
	UInt_t baseline_start = pretrigger - gatepre - sample_signal;

	TH1D *h_myoutput = new TH1D("h_myoutput", "Time difference", 1000, -20000, 20000);

	h_myoutput->SetFillColor(kBlue);

	TH1D *h_digi_FTS = new TH1D("h_digi_FTS", "Timestamp", 200, -1000, 3000);
	TH1D *h_my_FTS = new TH1D("h_my_FTS", "Timestamp", 200, -1000, 3000);
	h_digi_FTS->SetFillColor(kRed);
	h_my_FTS->SetFillColor(kBlue);

	auto loc_canv = new TCanvas("runs_from_delila", "c1", 1000, 800);
	loc_canv->Divide(3);
	Double_t sampleNumber_ch0[10000], rawSample_ch0[10000], sampleNumber_ch1[10000], rawSample_ch1[10000], timestamp[10000];
	std::vector<double> sampleNumber;
	entries = t_data->GetEntries();

	t_data->SetBranchStatus("Ch", kTRUE);
	t_data->SetBranchAddress("Ch", &t_channel);

	t_data->SetBranchStatus("RecordLength", kTRUE);
	t_data->SetBranchAddress("RecordLength", &recLength);

	t_data->SetBranchStatus("Signal", kTRUE);
	t_data->SetBranchAddress("Signal", data);

	double Tfine_ch0, Tfine_ch1, Ts0, Ts1;
	Long64_t digitmp_0, digitmp_1;

	for (j = 0; j < 10000; j++)
	{

		sampleNumber.push_back(j);
		timestamp[j] = j * 2;
	}

	Int_t sbzc0, sbzc1;
	double baseline_ch0;
	double baseline_ch1;
	Double_t cfd_delay0[10000], cfd_att0[10000], cfd0[10000], cfd_delay1[10000], cfd_att1[10000], cfd1[10000];

	TGraph *ch0, *pulse0;

	int iter = 0;
	for (int iEntry = 0; iEntry < 12001; iEntry++)
	{

		t_data->GetEntry(iEntry);

		if (t_channel == 2)
		{

			for (j = 0; j < recLength; j++)
			{

				sampleNumber_ch0[j] = j;
				rawSample_ch0[j] = data[j];
			}

			baseline_ch0 = 0;
			for (k = 0; k < 100; k++)
			{
				baseline_ch0 += rawSample_ch0[k];
			}
			baseline_ch0 /= 100;

			double cfd_delay0_new[10000] = {0};

			for (i = 0; i < recLength; i++)
			{
				cfd_att0[i] = ((rawSample_ch0[i] - baseline_ch0) * fraction_cfd);
				cfd_delay0[i] = -rawSample_ch0[i] + baseline_ch0;
			}
			for (i = delay_cfd; i < recLength; i++)
			{
				cfd_delay0_new[i + delay_cfd] = cfd_delay0[i];
			}
			for (i = delay_cfd; i < recLength; i++)
			{
				cfd0[i] = cfd_att0[i] + cfd_delay0_new[i];
			}

			Ts0 = Interpolation_orig(cfd0, sampleNumber, &Tfine_ch0);
		}

		if (t_channel == 3)
		{

			for (j = 0; j < recLength; j++)
			{
				sampleNumber_ch1[j] = j;
				rawSample_ch1[j] = data[j];
			}

			baseline_ch1 = 0;
			for (k = 0; k < 100; k++)
			{
				baseline_ch1 += rawSample_ch1[k];
			}
			baseline_ch1 /= 100;

			double cfd_delay1_new[10000] = {0};

			for (i = 0; i < recLength; i++)
			{
				cfd_att1[i] = ((rawSample_ch1[i] - baseline_ch1) * fraction_cfd);
				cfd_delay1[i] = -rawSample_ch1[i] + baseline_ch1;
			}
			for (i = delay_cfd; i < recLength; i++)
			{
				cfd_delay1_new[i + delay_cfd] = cfd_delay1[i];
			}
			for (i = delay_cfd; i < recLength; i++)
			{
				cfd1[i] = cfd_att1[i] + cfd_delay1_new[i];
			}

			Ts1 = Interpolation_orig(cfd1, sampleNumber, &Tfine_ch1);
		}

		if (t_channel == 2)
		{
			h_my_FTS->Fill(Tfine_ch0 * 2000);
		}

		h_myoutput->Fill(((Ts1 * 2000) - (Ts0 * 2000)));

		if (iter % 1000 == 0)
		{
			std::cout << "Ne aflam la evenimentul " << iter << endl;

			loc_canv->cd(1);
			h_myoutput->Draw();
			loc_canv->cd(2);
			h_my_FTS->Draw();
			loc_canv->cd(3);
			if (ch0)
				delete ch0;
			ch0 = new TGraph(recLength, timestamp, cfd0);
			ch0->SetTitle("CFD Signal");
			ch0->Draw("ALP");

			loc_canv->Update();
		}
		iter = iter + 1;
	}

	/*
		h_digi_output->Draw();
		h_myoutput->Draw("same");
		loc_canv->Update();
	*/

	/*
		loc_canv->cd(1);
		h_digi_output->Draw();
		loc_canv->cd(2);
		h_myoutput->Draw();

		loc_canv->cd(3);
		h_digi_FTS->Draw();
		loc_canv->cd(4);
		h_my_FTS->Draw();
		loc_canv->Update();
	*/
}