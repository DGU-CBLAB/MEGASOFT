////////////////////////////////////////////////////////////////////// 
// Metasoft.java
// (c) 2011-2020 Buhm Han
// 
// This file may be used for your personal use.
// This file may be modified but the modified version must retain this copyright notice.
// This file or modified version must not be redistributed
// without prior permission from the author.
// This software is provided “AS IS”, without warranty of any kind.
// In no event shall the author be liable for any claim, damages or other liability
// in connection with the software or the use or other dealings in the software.

import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;
import cern.jet.stat.Probability;

public class Metasoft
{
   // Arguments and default values
   private static String  inputFile_  = "";
   private static String  outputFile_ = "out";
   private static String  logFile_    = "log";
   private static String  pvalueTableFile_ = "HanEskinPvalueTable.txt";
   private static double  inputLambdaMeanEffect_    = 1.0;
   private static double  inputLambdaHeterogeneity_ = 1.0;
   private static boolean willComputeMvalue_ = false;
   private static double  priorSigma_ = 0.2;
   private static double  priorAlpha_ = 1.0;
   private static double  priorBeta_  = 1.0;
   private static double  mvaluePvalueThreshold_ = 1E-7;
   private static String  mvalueMethod_ = "exact";
   private static long    mcmcSample_ = 10000;
   private static long    mcmcBurnin_ = 1000;
   private static double  mcmcProbRandom_ = 0.01;
   private static double  mcmcMaxNumFlip_ = 0.1;
   private static boolean willComputeBinaryEffects_ = false;
   private static long    binaryEffectsSample_      = 1000;
   private static long    binaryEffectsLargeSample_ = 100000;
   private static double  binaryEffectsPvalueThreshold_ = 1E-4;
   private static int     seed_ = 0;
   private static boolean isVerbose_ = false;
   // Internally used variables
   private static int    numSnps_;
   private static int    maxNumStudy_;
   private static double outputLambdaMeanEffect_;
   private static double outputLambdaHeterogeneity_;
   private static ArrayList<Double> meanEffectParts_;
   private static ArrayList<Double> heterogeneityParts_;
   private static String argsSummary_;
   
   public static void main(String[] args) {
      double startTime = System.currentTimeMillis();
      handleArguments(args);
      System.err.printf("Arguments: "+argsSummary_);
      MetaSnp.readPvalueTableFile(pvalueTableFile_);
      System.err.println("----- Performing meta-analysis");
      doMetaAnalysis();
      computeLambda();
      printLog();
      System.err.println("----- Finished");
      double endTime   = System.currentTimeMillis();
      System.err.printf("----- Elapsed time: %.2f minutes\n", 
			(endTime - startTime)/(60*1000F));
   }
   
   private static void handleArguments(String[] args) {
      if (args.length == 0) {
	 printErrorAndQuit("ERROR: No argument. Please type 'java -jar Metasoft.jar -help' to see a list of options");
      }
      CommandLineParser parser = new GnuParser();
      Options options = new Options();
      // Option build-up
      options.addOption( OptionBuilder
			 .withLongOpt("input")
			 .withDescription("Input file (Required)")
			 .hasArg()
			 .withArgName("FILE")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("output")
			 .withDescription("Output file (default='out')")
			 .hasArg()
			 .withArgName("FILE")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("pvalue_table")
			 .withDescription("Pvalue table file (default='HanEskinPvalueTable.txt')")
			 .hasArg()
			 .withArgName("FILE")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("log")
			 .withDescription("Log file (default='log')")
			 .hasArg()
			 .withArgName("FILE")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("lambda_mean")
			 .withDescription("(Random Effects) User-specified lambda for mean effect part (default=1.0)")
			 .hasArg()
			 .withArgName("FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("lambda_hetero")
			 .withDescription("(Random Effects) User-specified lambda for heterogeneity part (default=1.0)")
			 .hasArg()
			 .withArgName("FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mvalue")
			 .withDescription("Compute m-value (default=false)")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mvalue_prior_sigma")
			 .withDescription("Sigma value for normal prior N(0, sigma^2) for effect (default=0.2)")
			 .hasArg()
			 .withArgName("FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mvalue_prior_beta")
			 .withDescription("Alpha and Beta value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
			 .hasArgs()
			 .withArgName("ALPHA BETA")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mvalue_p_thres")
			 .withDescription("Compute m-values only for SNPs whose FE or RE2 p-value is below this threshold (default=1E-7)")
			 .hasArg()
			 .withArgName("FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mvalue_method")
			 .withDescription("Which method to use to calculate m-value between 'exact' and 'mcmc' (default=exact)")
			 .hasArg()
			 .withArgName("METHOD_NAME")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mcmc_sample")
			 .withDescription("(MCMC) Number of samples (default=10,000)")
			 .hasArg()
			 .withArgName("INT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mcmc_burnin")
			 .withDescription("(MCMC) Number of burn-in (default=1,000)")
			 .hasArg()
			 .withArgName("INT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mcmc_prob_random_move")
			 .withDescription("(MCMC) Probability that a complete randomization move is suggested (default=0.01)")
			 .hasArg()
			 .withArgName("FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("mcmc_max_num_flip")
			 .withDescription("(MCMC) Usual move is flipping N bits where N ~ U(1,max_num_flip). If an integer value i >= 1 is given, max_num_flip = i. If a float value 0 < k < 1 is given, max_num_flip = k * #studies. (default=0.1)")
			 .hasArg()
			 .withArgName("INT or FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("binary_effects")
			 .withDescription("Compute binary effects model p-value (default=false)")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("binary_effects_sample")
			 .withDescription("(Binary effects) Number of importance sampling samples (default=1,000)")
			 .hasArg()
			 .withArgName("INT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("binary_effects_large")
			 .withDescription("(Binary effects) Large number of importance sampling samples for p-values above threshold (default=100,000)")
			 .hasArg()
			 .withArgName("INT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("binary_effects_p_thres")
			 .withDescription("(Binary effects) P-value threshold determining if we will use large number of samples (default=1E-4)")
			 .hasArg()
			 .withArgName("FLOAT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("seed")
			 .withDescription("Random number generator seed (default=0)")
			 .hasArg()
			 .withArgName("INT")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("verbose")
			 .withDescription("Print RSID verbosely per every 1,000 SNPs (default=false)")
			 .create() );
      options.addOption( OptionBuilder
			 .withLongOpt("help")
			 .withDescription("Print help")
			 .create() );
      // Parsing
      try {
	 CommandLine line = parser.parse(options, args);
	 inputFile_ = line.getOptionValue("input");
	 if (line.hasOption("output")) {
	    outputFile_ = line.getOptionValue("output");
	 }
	 if (line.hasOption("pvalue_table")) {
	    pvalueTableFile_ = line.getOptionValue("pvalue_table");
	 }
	 if (line.hasOption("log")) {
	    logFile_ = line.getOptionValue("log");
	 }
	 if (line.hasOption("lambda_mean")) {
	    inputLambdaMeanEffect_ = Double.valueOf(line.getOptionValue("lambda_mean"));
	 }
	 if (line.hasOption("lambda_hetero")) {
	    inputLambdaHeterogeneity_ = Double.valueOf(line.getOptionValue("lambda_hetero"));
	 }
	 if (line.hasOption("mvalue")) {
	    willComputeMvalue_ = true;
	    if (line.hasOption("mvalue_prior_sigma")) {
	       priorSigma_ = Double.valueOf(line.getOptionValue("mvalue_prior_sigma"));
	    }
	    if (line.hasOption("mvalue_prior_beta")) {
	       String[] c = line.getOptionValues("mvalue_prior_beta");
	       if (c.length < 2) {
		  printErrorAndQuit("Two arguments are needed for mvalue_prior_beta option");
	       }
	       priorAlpha_ = Double.valueOf(c[0]);
	       priorBeta_  = Double.valueOf(c[1]);
	    }
	    if (line.hasOption("mvalue_p_thres")) {
	       mvaluePvalueThreshold_ = Double.valueOf(line.getOptionValue("mvalue_p_thres"));
	    }
	    if (line.hasOption("mvalue_method")) {
	       mvalueMethod_ = line.getOptionValue("mvalue_method");
	    }
	    if (mvalueMethod_.equals("mcmc")) {
	       if (line.hasOption("mcmc_sample")) {
		  mcmcSample_ = Long.valueOf(line.getOptionValue("mcmc_sample"));
	       }
	       if (line.hasOption("mcmc_burnin")) {
		  mcmcBurnin_ = Long.valueOf(line.getOptionValue("mcmc_burnin"));
	       }
	       if (line.hasOption("mcmc_prob_random_move")) {
		  mcmcProbRandom_ = Double.valueOf(line.getOptionValue("mcmc_prob_random_move"));
	       }
	       if (line.hasOption("mcmc_max_num_flip")) {
		  mcmcMaxNumFlip_ = Double.valueOf(line.getOptionValue("mcmc_max_num_flip"));
	       }
	    }
	 }
	 if (line.hasOption("binary_effects")) {
	    willComputeBinaryEffects_ = true;
	    if (line.hasOption("binary_effects_sample")) {
	      binaryEffectsSample_ = Long.valueOf(line.getOptionValue("binary_effects_sample"));	       
	    }
	    if (line.hasOption("binary_effects_large")) {
	      binaryEffectsLargeSample_ = Long.valueOf(line.getOptionValue("binary_effects_large"));	       
	    }
	    if (line.hasOption("binary_effects_p_thres")) {
	      binaryEffectsPvalueThreshold_ = Double.valueOf(line.getOptionValue("binary_effects_p_thres"));	       
	    }
	 }
	 if (line.hasOption("seed")) {
	    seed_ = Integer.valueOf(line.getOptionValue("seed"));
	 }
	 if (line.hasOption("verbose")) {
	    isVerbose_ = true;
	 }
	 if (line.hasOption("help")) {
	    HelpFormatter formatter = new HelpFormatter();
	    formatter.setLongOptPrefix("-");
	    formatter.setWidth(100);
	    formatter.setOptionComparator(null);
	    formatter.printHelp("java -jar Metasoft.jar [options]",options);
	    System.out.println("------------------------------------------------");
	    System.out.println("The format of input_file:");
	    System.out.println("  Each row is each SNP.");
	    System.out.println("  1st column is RSID.");
	    System.out.println("  2nd and 3rd column are beta and its standard error for 1st study.");
	    System.out.println("  4th and 5th column are beta and its standard error for 2nd study.");
	    System.out.println("  6th and 7th column are beta and its standard error for 3rd study.");
	    System.out.println("  and so on...");
	    System.out.println("------------------------------------------------");
	    System.out.println();
	    System.exit(-1);
	 }
      } catch (ParseException exp) {
	 printErrorAndQuit(exp.getMessage()+
			   "\nPlease type 'java -jar Metasoft.jar -help' to see a list of options");
      }
      // Manual argument validity checkup
      System.err.println("-------- Processing arguments ---------");
      if (inputFile_.equals("")) {
	 printErrorAndQuit("A valid input file must be specified using option -input");
      }
      if (inputLambdaMeanEffect_ <= 0.0) {
	 printErrorAndQuit("lambda_mean option takes a float value > 0");
      }
      if (inputLambdaHeterogeneity_ <= 0.0) {
	 printErrorAndQuit("lambda_hetero option takes a float value > 0");
      }
      if (priorSigma_ <= 0.0) {
	 printErrorAndQuit("mvalue_prior_sigma option takes a float value > 0");
      }
      if (priorAlpha_ <= 0.0 || priorBeta_ <= 0.0) {
	 printErrorAndQuit("mvalue_prior_beta option takes two float values > 0");
      }
      if (mvaluePvalueThreshold_ < 0.0 || mvaluePvalueThreshold_ > 1.0) {
	 printErrorAndQuit("mvalue_p_thres takes a float value between 0 and 1");
      }
      if (!mvalueMethod_.equals("exact") &&
	  !mvalueMethod_.equals("mcmc")  &&
	  !mvalueMethod_.equals("variational")) {
	 printErrorAndQuit("mvalue_method option only takes a value 'exact' or 'mcmc'");
      }
      if (mcmcSample_ < 1) {
	 printErrorAndQuit("mcmc_sample option takes an integer value > 0");
      }
      if (mcmcBurnin_ < 1) {
	 printErrorAndQuit("mcmc_burnin option takes an integer value > 0");
      }
      if (mcmcSample_ < mcmcBurnin_) {
	 printErrorAndQuit("mcmc_sample must be larger than mcmc_burnin");
      }
      if (mcmcProbRandom_ < 0.0 || mcmcProbRandom_ > 1.0) {
	 printErrorAndQuit("mcmc_prob_random takes a float value between 0 and 1");
      }
      if (mcmcMaxNumFlip_ <= 0.0) {
	 printErrorAndQuit("mcmc_max_num_flip takes a value > 0");
      }
      if (binaryEffectsSample_ < 1) {
	 printErrorAndQuit("binary_effects_sample option takes an integer value > 0");
      }
      if (binaryEffectsLargeSample_ < 1) {
	 printErrorAndQuit("binary_effects_large option takes an integer value > 0");
      }
      if (binaryEffectsPvalueThreshold_ < 0.0 || binaryEffectsPvalueThreshold_ > 1.0) {
	 printErrorAndQuit("binary_effects_p_thres takes a float value between 0 and 1");
      }
      // Make summary for printing
      argsSummary_ = "";
      for (String s: args) {
	 argsSummary_ += " " + s;
      }
      argsSummary_ += "\n";
   }

   private static void printErrorAndQuit(String msg) {
      System.err.println(msg);
      System.exit(-1);
   }

   private static void doMetaAnalysis() {
      Random random = new Random(seed_);
      numSnps_ = 0;
      maxNumStudy_ = 0;
      meanEffectParts_    = new ArrayList<Double>();
      heterogeneityParts_ = new ArrayList<Double>();
      MetaSnp metaSnp;    // Store only 1 Snp at a time in memory.
      BufferedReader bufferedReader = null;
      PrintWriter printWriter = null;
      try {
	 bufferedReader = new BufferedReader(new FileReader(inputFile_));
      }
      catch (Exception exception) {
	 printErrorAndQuit("ERROR: Input file cannot be opened");
      }
      try {
	 printWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputFile_)));
      }
      catch (Exception exception) {
	 printErrorAndQuit("ERROR: Ouput file cannot be opened");
      }
      // Print headings
      MetaSnp.printHeadings(printWriter);
      try {
	 // Read 1 Snp information
	 String readLine;
	 while((readLine = bufferedReader.readLine()) != null) {
	    String[] tokens = readLine.split("\\s+");
	    if (tokens.length > 1) {             // only if non-empty
	       if (tokens[0].charAt(0) != '#') { // only if non-comment
		  String rsid = tokens[0];
		  metaSnp = new MetaSnp(rsid);
		  if (tokens.length % 2 == 0)
		     System.err.println("WARNING: # of Columns must be " +
					"odd including Rsid. Last column is ignored.");
		  int nStudy = (int)((tokens.length-1)/2);
		  if (nStudy > maxNumStudy_) {
		     maxNumStudy_ = nStudy;
		  }
		  for (int i = 0; i < nStudy; i++) {
		     double beta;
		     double standardError;
		     if (tokens[2*i+1].equals("NA") ||
			 tokens[2*i+1].equals("N/A") ||
			 tokens[2*i+2].equals("NA") ||
			 tokens[2*i+2].equals("N/A")) {
			metaSnp.addNaStudy();
		     } else {
			try {
			   beta = new Double(tokens[2*i+1]);
			   standardError = new Double(tokens[2*i+2]);
			   if (standardError <= 0.0) {
			      System.err.printf("Standard error cannot be <= zero "+
						"(%d th column is %f) "+
						"in the following line.\n",
						2*i+3, standardError);
			      System.err.println(readLine);
			      System.exit(-1);
			   }
			   metaSnp.addStudy(beta, standardError);
			} 
			catch (Exception exception) {
			   System.err.println("Incorrect float value in following line. "+
					      "Possibly not a double");
			   System.err.println(readLine);
			   System.exit(-1);
			}
		     }
		  }
		  if (metaSnp.getNStudy() > 1) {
		     // Analyze 1 Snp on-the-fly.
		     if (isVerbose_ && numSnps_ % 1000 == 0) {
			System.err.printf("Analyzing SNP #%d (%s)\n", numSnps_ + 1, rsid);
		     }
		     // FE, RE, and New RE
		     metaSnp.computeFixedEffects();
		     metaSnp.computeRandomEffects();
		     metaSnp.computeHanEskin(inputLambdaMeanEffect_,
					     inputLambdaHeterogeneity_);
		     meanEffectParts_.add(metaSnp.getStatisticHanEskinMeanEffectPart());
		     double h = metaSnp.getStatisticHanEskinHeterogeneityPart();
		     if (h > 0.0) {
			heterogeneityParts_.add(h);
		     }
		     // Binary effects model
		     if (willComputeBinaryEffects_) {
			metaSnp.computeBinaryEffectsPvalue(binaryEffectsSample_, random.nextInt());
			if (metaSnp.getPvalueBinaryEffects() <= binaryEffectsPvalueThreshold_) {
			   metaSnp.computeBinaryEffectsPvalue(binaryEffectsLargeSample_, random.nextInt());
			}
		     }
		     // Mvalues
		     if (willComputeMvalue_) {
			if (metaSnp.getPvalueFixedEffects() <= mvaluePvalueThreshold_ ||
			    metaSnp.getPvalueHanEskin()     <= mvaluePvalueThreshold_) {
			   if (mvalueMethod_.equals("exact")) {
			      metaSnp.computeMvalues(priorAlpha_, priorBeta_, priorSigma_);
			   } else if (mvalueMethod_.equals("mcmc")) {
			      metaSnp.computeMvaluesMCMC(priorAlpha_, priorBeta_, priorSigma_,
							 mcmcSample_, mcmcBurnin_, mcmcProbRandom_, 
							 mcmcMaxNumFlip_, 
							 random.nextInt());
			   } else {
			      assert false : mvalueMethod_;
			   }
			}
		     }
		     numSnps_++;
		  }
		  metaSnp.printResults(printWriter);

	       }
	    }
	 }
      }
      catch (IOException exception) {
	 System.err.println("ERROR: error encountered while reading input file");
	 System.exit(-1);
      }
      try {
	 bufferedReader.close();
	 printWriter.close();
      }
      catch (Exception exception) {
	 System.err.println("ERROR: file cannot be closed");
	 System.exit(-1);
      }
   }
   
   private static void computeLambda() {
      double median;
      double expectedMedian;
      if (meanEffectParts_.size() > 0) {
	 Collections.sort(meanEffectParts_);
	 median = meanEffectParts_.get((int)(meanEffectParts_.size()/2.0));
	 expectedMedian = Math.pow(Probability.normalInverse(0.25),2.0);
	 outputLambdaMeanEffect_ = median / expectedMedian;
      }
      if (heterogeneityParts_.size() > 0) {
	 Collections.sort(heterogeneityParts_);
	 median = heterogeneityParts_.get((int)(heterogeneityParts_.size()/2.0));
	 if (maxNumStudy_ > 50) {
	    expectedMedian = Math.pow(Probability.normalInverse(0.25),2.0);
	 } else {
	    expectedMedian = expectedMedianHanEskinHeterogeneityPart_[maxNumStudy_ - 2];
	 }
	 outputLambdaHeterogeneity_ = median / expectedMedian;
      }
   }

   private static void printLog() {
      try {
	 PrintWriter printWriter = new PrintWriter(new BufferedWriter(new FileWriter(logFile_)));
	 printWriter.printf("Arguments: %s ", argsSummary_);
	 printWriter.printf("Input File: %s\n", inputFile_);
	 printWriter.printf("Output File: %s\n", outputFile_);
	 printWriter.printf("Log File: %s\n", logFile_);
	 printWriter.printf("p-value Table File: %s\n", pvalueTableFile_);
	 printWriter.printf("Number of SNPs analyzed: %d\n", numSnps_);
	 printWriter.printf("Maximum number of studies: %d\n", maxNumStudy_);
	 printWriter.printf("Specified lambda for   mean effect part (default = 1.0): %f\n", inputLambdaMeanEffect_);
	 printWriter.printf("Specified lambda for heterogeneity part (default = 1.0): %f\n", inputLambdaHeterogeneity_);
	 printWriter.printf("Newly calculated inflation factor lambda for   mean effect part: %f\n", outputLambdaMeanEffect_);
	 printWriter.printf("Newly calculated inflation factor lambda for heterogeneity part: %f\n", outputLambdaHeterogeneity_);
	 printWriter.close();
      }
      catch (Exception exception) {
	 System.err.println("ERROR: error encountered while writing in log file");
	 System.exit(-1);
      }
   }

   private static double[] expectedMedianHanEskinHeterogeneityPart_ = // from nStudy 2 to 50
   {0.2195907137,0.2471516439,0.2642270318,0.2780769264,0.2886280267,0.2977812664,0.3020051148,0.3091428179,0.3158605559,0.3221788173,0.3259133140,0.3295976587,0.3335375196,0.3358395088,0.3368309971,0.3421941686,0.3448030927,0.3463590948,0.3477384754,0.3487900288,0.3494938171,0.3542087791,0.3573286353,0.3589703411,0.3586951356,0.3596101209,0.3605611682,0.3624799993,0.3648322669,0.3659817739,0.3671267389,0.3693952373,0.3693395144,0.3696863113,0.3706067524,0.3718103285,0.3749536619,0.3758886239,0.3753612342,0.3781458299,0.3798346038,0.3763434983,0.3796968747,0.3784334922,0.3794411347,0.3808582942,0.3813485882,0.3843230993,0.3824863479};
}
