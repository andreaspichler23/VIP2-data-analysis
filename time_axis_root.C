{
  TCanvas* c1 = new TCanvas("c1", "A graph", 200, 10, 700, 500);

  // make up some data
  const int n = 20;
  char label[n][50];
  double y[n];

  for (int i=0; i<n; i++) {
    float xx = i * 0.1;
    y[i] = 10 * sin(xx + 0.2);
    sprintf(label[i], "2012-01-17 09:%02i:00", i);
  }

  // convert (t,y) into (x,y) for drawing
  double x[n];
  for (int i=0; i<n; i++) {
    TDatime d(label[i]); // format of label must be "yyyy-mm-dd hh:mm:ss" here
    x[i] = d.Convert();  // convert to seconds
    cout << label[i] << " " << x[i] << endl;
  }

  // draw a graph with (x,y) coords
  TGraph* gr = new TGraph(n, x, y);
  gr->GetXaxis()->SetTitle("Time");
  gr->GetYaxis()->SetTitle("Sunspots");
  gr->Draw("AL*");

  // set x axis to display time
  gr->GetXaxis()->SetTimeDisplay(1);

  // set the format you want displayed
  gr->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M:%S"); // normal
  // gr->GetXaxis()->SetTimeFormat("#splitline{%Y-%m-%d}{%H:%M:%S}");    // fancy

  // provide an offset, added to the entire series
  gr->GetXaxis()->SetTimeOffset(0);

  // tweak label appearance
  //gr->GetXaxis()->SetLabelSize(0.05);

  // draw!
  gr->Draw("AL*");

  return c1;
}
