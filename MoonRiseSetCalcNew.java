import java.util.Date;


public class MoonRiseSetCalcNew {

	double  PI = Math.PI;
	double DR = PI/180;
	double K1 = 15*DR*1.0027379;

	boolean Moonrise = false;
	boolean Moonset  = false;

	int[] Rise_time = {0, 0};
	int[] Set_time  = {0, 0};
	double Rise_az = 0.0;
	double Set_az  = 0.0;

	double[] Sky = {0.0, 0.0, 0.0};
	double[] RAn = {0.0, 0.0, 0.0};
	double[] Dec = {0.0, 0.0, 0.0};
	double[] VHz = {0.0, 0.0, 0.0};

	Date Now ;

	String thisday;
	String moonrise= "";
	String moonset = "";
	String azRise = "";
	String azSet = "";
	
	

	
	double timeZone ;
	
	boolean makeTimeZoneSignCurrection;
	 
	public MoonRiseSetCalcNew(double userTimeZone, double userDst){
		
		
		
		timeZone =  -(userTimeZone)*1.0 + userDst;
		
		if(userTimeZone>=0){
			makeTimeZoneSignCurrection = true;
		}else{
			makeTimeZoneSignCurrection = false;
		}
		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		// user input for New Delhi for 2 Nov 2012
		
		double userTimeZone = 5.5, userDst = 0.0;
		
		int lat_deg = 28;
		int lat_min = 36;
		int lat_sec = 0;
		boolean south = false; 
		
		int lon_deg = 77;  
		int lon_min = 12;
		int lon_sec = 0;
		boolean west = false;
		
		int day = 2;
		int month =11;
		int year = 2012;
		
		new MoonRiseSetCalcNew(userTimeZone, userDst).compute(lat_deg,lat_min ,lat_sec ,south ,lon_deg ,lon_min ,lon_sec ,west, day,month,year);
	}

	private void compute(int lat_degrees,int lat_minutes,int lat_secundes,boolean south,int lon_degrees,int lon_minutes,int lon_secundes,boolean west, int day, int month, int year)
	{
				
		
		month = month-1 ;
		
		 double lat = lat_degrees + lat_minutes/60.0 + lat_secundes/3600.0;
		 double lon = lon_degrees + lon_minutes/60.0 + lon_secundes/3600.0;

	    if (south) lat = -lat;
	    if (west) lon = -lon;

	    Now = new Date(year, month, day, 0, 0, 0);
	    showdate(Now);

	    moonrise = "";
	    moonset  = "";
	    riseset(lat, lon);

	    showResult();
	}
	private void showResult() {
		System.out.println("Moon Rise : "+moonrise);
		System.out.println("Moon Set : "+moonset);
		System.out.println("azRise : "+azRise);
		System.out.println("azSet : "+azSet);
		System.out.println("Date : "+thisday);
	}
	// display the date
	private void showdate(Date d)
	{
	    int date = d.getDate(); 
	    String thisday   = d.getDate() + " ";
	    int thismonth = d.getMonth()+1;
	    String thisyear  = " " + d.getYear();
		thisday = thisday + thismonth + thisyear;
	    this.thisday = thisday;
	}
	// calculate moonrise and moonset times
	private void riseset(double lat,double lon )
	{
	    int i, j, k;

	    double jdlp = julian_day();			// stored for Lunar Phase calculation
	    double jd = jdlp - 2451545;           // Julian day relative to Jan 1.5, 2000
	    double[][] mp = new double[3][3];                     // create a 3x3 array
	    for (i = 0; i < 3; i++)
	    {
	      for (j = 0; j < 3; j++)
	            mp[i][j] = 0.0;
	    }

	    lon = lon/360;

	    double tz =Math.abs(timeZone/24); 
	    
	    if(makeTimeZoneSignCurrection)
	    	tz = -tz;
	    double t0 = lst(lon, jd, tz);                 // local sidereal time

	    jd = jd + tz;                              // get moon position at start of day

	    for (k = 0; k < 3; k++)
	    {   
	    	moon(jd);
	        mp[k][0] = Sky[0];
	        mp[k][1] = Sky[1];
	        mp[k][2] = Sky[2];
	        jd = jd + 0.5;      
	    }   

	    if (mp[1][0] <= mp[0][0])
	        mp[1][0] = mp[1][0] + 2*PI;

	    if (mp[2][0] <= mp[1][0])
	        mp[2][0] = mp[2][0] + 2*PI;

	    RAn[0] = mp[0][0];
	    Dec[0] = mp[0][1];

	    Moonrise = false;                          // initialize
	    Moonset  = false;
	    double ph = 0.0;
	    double temp;
	    for (k = 0; k < 24; k++)                   // check each hour of this day
	    {
	    	temp = k;
	    	ph = (temp + 1)/24;
	        
	        RAn[2] = interpolate(mp[0][0], mp[1][0], mp[2][0], ph);
	        Dec[2] = interpolate(mp[0][1], mp[1][1], mp[2][1], ph);
	        
	        VHz[2] = test_moon(k, timeZone, t0, lat, mp[1][2]);

	        RAn[0] = RAn[2];                       // advance to next hour
	        Dec[0] = Dec[2];
	        VHz[0] = VHz[2];
	    }

	    // display results

	    moonrise = zintstr(Rise_time[0], 2) + ":" + zintstr(Rise_time[1], 2);
	    moonset  = zintstr( Set_time[0], 2) + ":" + zintstr( Set_time[1], 2);


	        
	if (Moonrise) {azRise = frealstr(Rise_az, 5, 1) + "°";} else {azRise = " ";}
	if (Moonset)  {azSet  = frealstr(Set_az, 5, 1) + "°";}  else {azSet = " ";}

	    special_message();
	}

	// check for no moonrise and/or no moonset
	private void special_message()
	{
	    if ((!Moonrise)&&(!Moonset))               // neither moonrise nor moonset
	    {
	        if (VHz[2] < 0)
	           moonrise = "Moon down all day";
	        else
	            moonrise = "Moon up all day";

	        moonset = "";
	    }
	    else                                       // moonrise or moonset
	    {
	        if (!Moonrise)
	            moonrise = "No moonrise this date";
	        else if (!Moonset)
	           moonset  = "No moonset this date";
	    }
	}
	private double julian_day()							// be carefull, the function of the similare name (Julian_Day) is used in astroAWK1.js library
	// current function uses date in the form of "date object", while that other function uses three arguments of calendar date
{
		double a, b, jd;
		boolean gregorian;
		
		int month = Now.getMonth() + 1;
		int day   = Now.getDate();
		int year  = Now.getYear();
		
		gregorian = (year < 1583) ? false : true;
		
		if ((month == 1)||(month == 2))
		{
		year  = year  - 1;
		month = month + 12;
		}
		
		a = Math.floor(year/100);
		if (gregorian) b = 2 - a + Math.floor(a/4);
		else           b = 0.0;
		
		jd = Math.floor(365.25*(year + 4716)) 
		+ Math.floor(30.6001*(month + 1)) 
		+ day + b - 1524.5;
		
		return jd;
}
	// Local Sidereal Time for timeZone
	private double lst( double lon, double jd, double z )
	{
		double s = 24110.5 + 8640184.812999999*jd/36525 + 86636.6*z + 86400*lon;
	    s = s/86400;
	    s = s - Math.floor(s);
	    return s*360*DR;
	}
	
	// moon's position using fundamental arguments 
	// (Van Flandern & Pulkkinen, 1979)
	private void moon( double jd )
	{
	    double d, f, g, h, m, n, s, u, v, w;

	    h = 0.606434 + 0.03660110129*jd;
	    m = 0.374897 + 0.03629164709*jd;
	    f = 0.259091 + 0.0367481952 *jd;
	    d = 0.827362 + 0.03386319198*jd;
	    n = 0.347343 - 0.00014709391*jd;
	    g = 0.993126 + 0.0027377785 *jd;

	    h = h - Math.floor(h);
	    m = m - Math.floor(m);
	    f = f - Math.floor(f);
	    d = d - Math.floor(d);
	    n = n - Math.floor(n);
	    g = g - Math.floor(g);

	    h = h*2*PI;
	    m = m*2*PI;
	    f = f*2*PI;
	    d = d*2*PI;
	    n = n*2*PI;
	    g = g*2*PI;

	    v = 0.39558*Math.sin(f + n);
	    v = v + 0.082  *Math.sin(f);
	    v = v + 0.03257*Math.sin(m - f - n);
	    v = v + 0.01092*Math.sin(m + f + n);
	    v = v + 0.00666*Math.sin(m - f);
	    v = v - 0.00644*Math.sin(m + f - 2*d + n);
	    v = v - 0.00331*Math.sin(f - 2*d + n);
	    v = v - 0.00304*Math.sin(f - 2*d);
	    v = v - 0.0024 *Math.sin(m - f - 2*d - n);
	    v = v + 0.00226*Math.sin(m + f);
	    v = v - 0.00108*Math.sin(m + f - 2*d);
	    v = v - 0.00079*Math.sin(f - n);
	    v = v + 0.00078*Math.sin(f + 2*d + n);
	    
	    u = 1 - 0.10828*Math.cos(m);
	    u = u - 0.0188 *Math.cos(m - 2*d);
	    u = u - 0.01479*Math.cos(2*d);
	    u = u + 0.00181*Math.cos(2*m - 2*d);
	    u = u - 0.00147*Math.cos(2*m);
	    u = u - 0.00105*Math.cos(2*d - g);
	    u = u - 0.00075*Math.cos(m - 2*d + g);
	    
	    w = 0.10478*Math.sin(m);
	    w = w - 0.04105*Math.sin(2*f + 2*n);
	    w = w - 0.0213 *Math.sin(m - 2*d);
	    w = w - 0.01779*Math.sin(2*f + n);
	    w = w + 0.01774*Math.sin(n);
	    w = w + 0.00987*Math.sin(2*d);
	    w = w - 0.00338*Math.sin(m - 2*f - 2*n);
	    w = w - 0.00309*Math.sin(g);
	    w = w - 0.0019 *Math.sin(2*f);
	    w = w - 0.00144*Math.sin(m + n);
	    w = w - 0.00144*Math.sin(m - 2*f - n);
	    w = w - 0.00113*Math.sin(m + 2*f + 2*n);
	    w = w - 0.00094*Math.sin(m - 2*d + g);
	    w = w - 0.00092*Math.sin(2*m - 2*d);

	    s = w/Math.sqrt(u - v*v);                  // compute moon's right ascension ...  
	    Sky[0] = h + Math.atan(s/Math.sqrt(1 - s*s));

	    s = v/Math.sqrt(u);                        // declination ...
	    Sky[1] = Math.atan(s/Math.sqrt(1 - s*s));

	    Sky[2] = 60.40974*Math.sqrt( u );          // and parallax
	}

	// 3-point interpolation
	private double interpolate( double f0, double f1, double f2,double p )
	{
		double a = f1 - f0;
		double b = f2 - f1 - a;
		double f = f0 + p*(2*a + b*(2*p - 1));

	    return f;
	}
	
	// test an hour for an event
	private double test_moon( int k, double zone, double t0, double lat, double plx )
	{
	    double[] ha = {0.0, 0.0, 0.0};
	    double a, b, c, d, e, s, z;
	    int hr, min;
	    double time;
	    double az, hz, nz, dz;

	    if (RAn[2] < RAn[0])
	        RAn[2] = RAn[2] + 2*PI;
	    
	    ha[0] = t0 - RAn[0] + k*K1;
	    ha[2] = t0 - RAn[2] + k*K1 + K1;
	    
	    ha[1]  = (ha[2] + ha[0])/2;                // hour angle at half hour
	    Dec[1] = (Dec[2] + Dec[0])/2;              // declination at half hour

	    s = Math.sin(DR*lat);
	    c = Math.cos(DR*lat);

	    // refraction + sun semidiameter at horizon + parallax correction
	    z = Math.cos(DR*(90.567 - 41.685/plx));

	    if (k <= 0)                                // first call of function
	        VHz[0] = s*Math.sin(Dec[0]) + c*Math.cos(Dec[0])*Math.cos(ha[0]) - z;

	    VHz[2] = s*Math.sin(Dec[2]) + c*Math.cos(Dec[2])*Math.cos(ha[2]) - z;
	    
	    if (sgn(VHz[0]) == sgn(VHz[2]))
	        return VHz[2];                         // no event this hour
	    
	    VHz[1] = s*Math.sin(Dec[1]) + c*Math.cos(Dec[1])*Math.cos(ha[1]) - z;

	    a = 2*VHz[2] - 4*VHz[1] + 2*VHz[0];
	    b = 4*VHz[1] - 3*VHz[0] - VHz[2];
	    d = b*b - 4*a*VHz[0];

	    if (d < 0)
	        return VHz[2];                         // no event this hour
	    
	    d = Math.sqrt(d);
	    e = (-b + d)/(2*a);

	    if (( e > 1 )||( e < 0 ))
	        e = (-b - d)/(2*a);

	    time = k + e + 1/120;                      // time of an event + round up
	    hr   = (int) Math.floor(time);
	    min  = (int) Math.floor((time - hr)*60);

	    hz = ha[0] + e*(ha[2] - ha[0]);            // azimuth of the moon at the event
	    nz = -Math.cos(Dec[1])*Math.sin(hz);
	    dz = c*Math.sin(Dec[1]) - s*Math.cos(Dec[1])*Math.cos(hz);
	    az = Math.atan2(nz, dz)/DR;
	    if (az < 0) az = az + 360;
	    
	    if ((VHz[0] < 0)&&(VHz[2] > 0))
	    {
	        Rise_time[0] = hr;
	        Rise_time[1] = min;
	        Rise_az = az;
	        Moonrise = true;
	    }
	    
	    if ((VHz[0] > 0)&&(VHz[2] < 0))
	    {
	        Set_time[0] = hr;
	        Set_time[1] = min;
	        Set_az = az;
	        Moonset = true;
	    }

	    return VHz[2];
	}
	// returns value for sign of argument
	private int sgn(double x )
	{
	    int  rv;
	    if (x > 0.0)      rv =  1;
	    else if (x < 0.0) rv = -1;
	    else              rv =  0;
	    return rv;
	}
	// format a positive integer with leading zeroes
	private String zintstr( int num, int width )
	{
//	    String str = num.toString(10);
		String str = String.valueOf(num);
	    int len = str.length();
	    String intgr = "";
	    int i;

	    for (i = 0; i < width - len; i++)          // append leading zeroes
	        intgr += '0';

	    for (i = 0; i < len; i++)                  // append digits
	        intgr += str.charAt(i);

	    return intgr;
	}
	// format a real number
	private String frealstr( double num, int width, int fract )
	{
	    /*var str = num.toFixed(fract);
	    var len = str.length;
	    var real = "";
	    var i;*/
		num = Math.round( num * 100.0 ) / 100.0;
		String str = String.valueOf(num);
//	    var str = num.toFixed(fract);
	    int len = str.length();
	    String real = "";
	    int i;

	    for (i = 0; i < width - len; i++)          // append leading spaces
	        real += ' ';

	    for (i = 0; i < len; i++)                  // append digits
	        real += str.charAt(i);

	    return real;
	}
}
