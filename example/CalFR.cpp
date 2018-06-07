
double Pow2(double x) { return TMath::Power(x, 2.0); }

void FrHistoCal(TH1F* hfrac1, TH1F* hfrac2, TH1F* hfr1, TH1F* hfr2, TH1F* hfb, TH1F* hfl){

  /*
   calculate the b jet and light jet fakerate using sample 1 and 2, 
    where hfrac1, hfrac2 is b jet fraction distribution in sample1, 2
          hfr1, hfr2 is fakerate distribution in sample1, 2
    results saved in hfb and hfl
  */

  for(int i=1; i<= hfr1->GetNbinsX(); i++){

    int ibinfrac = hfrac1->FindBin(hfr1->GetBinCenter(i));

    double fb1 = hfrac1->GetBinContent(ibinfrac);           double err_fb1 = hfrac1->GetBinError(ibinfrac);
    double fl1 = 1.0 - fb1;
    double fb2 = hfrac2->GetBinContent(ibinfrac);           double err_fb2 = hfrac2->GetBinError(ibinfrac);
    double fl2 = 1.0 - fb2;
    double FR1 = hfr1->GetBinContent(i);                    double err_FR1 = hfr1->GetBinError(i);
    double FR2 = hfr2->GetBinContent(i);                    double err_FR2 = hfr2->GetBinError(i);

    double norm = 1.0/(fb1- fb2);

    double FR_b = norm*( fl2*FR1 - fl1*FR2); // fakerate for b jets
    double FR_l = norm*(-fb2*FR1 + fb1*FR2); // fakerate for light jets

    double coef_b_fb1 = -fl2*(FR2-FR1)*Pow2(norm);
    double coef_b_fb2 =  fl1*(FR2-FR1)*Pow2(norm);
    double coef_b_FR1 = -fl2*norm;
    double coef_b_FR2 =  fl1*norm;

    double coef_l_fb1 = -fb2*(FR2-FR1)*Pow2(norm);
    double coef_l_fb2 =  fb1*(FR2-FR1)*Pow2(norm);
    double coef_l_FR1 = -fb2*norm;
    double coef_l_FR2 =  fb1*norm;

    // error for fakerate of b jets
    double err_FR_b  = ( norm*( fl2*FR1 - fl1*FR2)>=0.0 ?  TMath::Sqrt(Pow2(coef_b_fb1*err_fb1) + Pow2(coef_b_fb2*err_fb2) + Pow2(coef_b_FR1*err_FR1) + Pow2(coef_b_FR2*err_FR2)) : 0.0 );
    // error for fakerate of light jets
    double err_FR_l  = ( norm*(-fb2*FR1 + fb1*FR2)>=0.0 ?  TMath::Sqrt(Pow2(coef_l_fb1*err_fb1) + Pow2(coef_l_fb2*err_fb2) + Pow2(coef_l_FR1*err_FR1) + Pow2(coef_l_FR2*err_FR2)) : 0.0 );

    FR_b = ( FR_b >= 0.0 ? FR_b : 0.0 );
    FR_l = ( FR_l >= 0.0 ? FR_l : 0.0 );

    hfb->SetBinContent( i, FR_b);   hfb->SetBinError( i, err_FR_b);
    hfl->SetBinContent( i, FR_l);   hfl->SetBinError( i, err_FR_l);

    std::cout << " Bin Number " << i << " calc B " << std::setw(10) << hfb->GetBinContent(i) << " err " << std::setw(10) << hfb->GetBinError(i) << std::endl;
    std::cout << " Bin Number " << i << " calc L " << std::setw(10) << hfl->GetBinContent(i) << " err " << std::setw(10) << hfl->GetBinError(i) << std::endl;
  }

}

void FrHistoCalTot(TH1F* hfrac1, TH1F* hfrac2, TH1F* hfr1, TH1F* hfr2, double fbtot, double err_fbtot, TH1F* hftot){

  /*
    calculate the averaged fakerate hftot
    where hfrac1, hfrac2 is b jet fraction distribution in GJet sample1, 2
          hfr1, hfr2 is fakerate distribution in GJet sample1, 2
          fbtot, err_fbtot is the b jet fraction and error in the 0tag/1tag QCD sample
    results saved in hftot

    hftot will be used for closure test or background prediction
  */

  for(int i=1; i<= hfr1->GetNbinsX(); i++){

    int ibinfrac = hfrac1->FindBin(hfr1->GetBinCenter(i));

    double fb1 = hfrac1->GetBinContent(ibinfrac);           double err_fb1 = hfrac1->GetBinError(ibinfrac);
    double fb2 = hfrac2->GetBinContent(ibinfrac);           double err_fb2 = hfrac2->GetBinError(ibinfrac);
    double FR1 = hfr1->GetBinContent(i);                    double err_FR1 = hfr1->GetBinError(i);
    double FR2 = hfr2->GetBinContent(i);                    double err_FR2 = hfr2->GetBinError(i);

    double norm = 1.0/(fb1- fb2);

    // averaged fakerate
    double FRtot = norm*( (fbtot-fb2)*FR1 + (fb1-fbtot)*FR2 );
    
    double coef_FR1   = (fbtot - fb2)*norm; 
    double coef_FR2   = (fb1 - fbtot)*norm;
    double coef_fbtot = (FR1 - FR2)  *norm;
    double coef_fb1   = (fbtot - fb2)*(FR2-FR1)*norm*norm;
    double coef_fb2   = (fbtot - fb1)*(FR1-FR2)*norm*norm;

    // uncertainty of the averaged fakerate
    double err_FRtot = TMath::Sqrt(Pow2(coef_FR1*err_FR1) + Pow2(coef_FR2*err_FR2) + Pow2(coef_fbtot*err_fbtot) + Pow2(coef_fb1*err_fb1) + Pow2(coef_fb2*err_fb2));

    hftot->SetBinContent(i, FRtot); hftot->SetBinError(i, err_FRtot);
   
    std::cout << " Bin Number " << i << " calc fakerate " << std::setw(10) << hftot->GetBinContent(i) << " err " << std::setw(10) << hftot->GetBinError(i) << std::endl;
  }
}
