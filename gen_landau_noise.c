

void gen_landau_noise( int ngen=30000000, float mu=0., float sigma=0.05, int seed=12345 ) {


	char outfilename[1000] ;
	sprintf( outfilename, "landau-noise--mu%.2f--sigma%.2f--seed%d.root", mu, sigma, seed ) ;
	printf("\n\n Output file: %s\n\n", outfilename ) ;

	TRandom2 ran = TRandom2( seed ) ;

	TFile* tf_out = new TFile( outfilename, "recreate" ) ;
	char title[1000] ;
	sprintf( title, "Landau, mu=%.4f, sigma=%.4f, seed=%d", mu, sigma, seed ) ;
        TTree* tt_out = new TTree( "rantree", title ) ;

	float dpx, dpy, dpz ;
	tt_out -> Branch( "dpx", &dpx, "dpx/F" ) ;
	tt_out -> Branch( "dpy", &dpy, "dpy/F" ) ;
	tt_out -> Branch( "dpz", &dpz, "dpz/F" ) ;

	for ( int i=0; i<ngen; i++ ) {
		if ( i%100 == 0 ) printf("   %9d / %9d\r", i, ngen ) ;
		dpx = ran.Landau( mu, sigma ) ;
		dpy = ran.Landau( mu, sigma ) ;
		dpz = ran.Landau( mu, sigma ) ;
		if ( ran.Uniform() < 0.5 ) dpx = -1. * dpx ;
		if ( ran.Uniform() < 0.5 ) dpy = -1. * dpy ;
		if ( ran.Uniform() < 0.5 ) dpz = -1. * dpz ;
		tt_out -> Fill() ;
	} // i

	printf("\n\n Done.\n\n") ;

	tt_out->Write() ;

	tf_out->Close() ;



}

