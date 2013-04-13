// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Dirichlet.SampleFromPseudoCounts"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "sample", "pseudoCounts" }, typeof(Dirichlet), "SampleFromPseudoCounts")]
	[Quality(QualityBand.Stable)]
	public static class DirichletFromPseudoCountsOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sampleFromPseudoCounts'.</param>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sampleFromPseudoCounts) p(sampleFromPseudoCounts) factor(sampleFromPseudoCounts,pseudoCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Dirichlet sample, Vector pseudoCounts, [Fresh] Dirichlet to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sampleFromPseudoCounts'.</param>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sampleFromPseudoCounts) p(sampleFromPseudoCounts) factor(sampleFromPseudoCounts,pseudoCount) / sum_sampleFromPseudoCounts p(sampleFromPseudoCounts) messageTo(sampleFromPseudoCounts))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Dirichlet sample, Vector pseudoCounts) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sampleFromPseudoCounts'.</param>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sampleFromPseudoCounts) p(sampleFromPseudoCounts) log(factor(sampleFromPseudoCounts,pseudoCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Dirichlet sample, Vector pseudoCounts, [Fresh] Dirichlet to_sample)
		{
			return to_sample.GetAverageLog(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sampleFromPseudoCounts'.</param>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sampleFromPseudoCounts,pseudoCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector sample, Vector pseudoCounts)
		{
			return SampleAverageConditional(pseudoCounts).GetLogProb(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sampleFromPseudoCounts'.</param>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sampleFromPseudoCounts,pseudoCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector sample, Vector pseudoCounts)
		{
			return LogAverageFactor(sample, pseudoCounts);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sampleFromPseudoCounts'.</param>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sampleFromPseudoCounts,pseudoCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Vector sample, Vector pseudoCounts)
		{
			return LogAverageFactor(sample, pseudoCounts);
		}

		/// <summary>
		/// EP message to 'sampleFromPseudoCounts'
		/// </summary>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <returns>The outgoing EP message to the 'sampleFromPseudoCounts' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sampleFromPseudoCounts' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet SampleAverageConditional(Vector pseudoCounts)
		{
			return new Dirichlet(pseudoCounts);
		}
		/// <summary>
		/// VMP message to 'sampleFromPseudoCounts'
		/// </summary>
		/// <param name="pseudoCounts">Constant value for 'pseudoCount'.</param>
		/// <returns>The outgoing VMP message to the 'sampleFromPseudoCounts' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sampleFromPseudoCounts' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet SampleAverageLogarithm(Vector pseudoCounts)
		{
			return new Dirichlet(pseudoCounts);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.DirichletSymmetric"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "DirichletSymmetric")]
	[Quality(QualityBand.Preview)]
	[Buffers("probMeanLog")]
	public static class DirichletSymmetricOp
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="alpha">Incoming message from 'alpha'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="probMeanLog">Buffer for E[log(prob)]</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(prob,alpha) p(prob,alpha) log(factor(prob,alpha,K))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="alpha"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static double AverageLogFactor([Proper] Dirichlet prob, [SkipIfUniform] Gamma alpha, [Fresh] Vector probMeanLog)
		{
			if (alpha.IsPointMass)
				return AverageLogFactor(prob, alpha.Point, probMeanLog);
			double SumElogP = probMeanLog.Sum();
			double K = (double)probMeanLog.Count;
			double a = alpha.Shape;
			double b = alpha.Rate;
			double averageFactor = GammaFromShapeAndRateOp.ELogGamma(Gamma.FromShapeAndRate(a, b / K));
			averageFactor -= K * GammaFromShapeAndRateOp.ELogGamma(Gamma.FromShapeAndRate(a, b));
			averageFactor += (a / b - 1.0) * SumElogP;
			return averageFactor;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="alpha">Incoming message for 'alpha'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>int p(alpha) log(factor(prob,alpha,K))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Quality(QualityBand.Experimental)]
		public static double AverageLogFactor(Vector prob, [SkipIfUniform] Gamma alpha)
		{
			Dirichlet d = Dirichlet.PointMass(prob);
			return AverageLogFactor(d, alpha, d.GetMeanLog());
		}

        /// <summary>
        /// Evidence message for VMP
        /// </summary>
        /// <param name="prob">Constant value for 'prob'.</param>
        /// <param name="alpha">Incoming message for 'alpha'.</param>
        /// <returns>Average of the factor's log-value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>int p(alpha) log(factor(prob,alpha,K))</c>.
        /// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
        /// </para></remarks>
        [Quality(QualityBand.Experimental)]
        public static double AverageLogFactor(Vector prob, [SkipIfUniform] ConjugateDirichlet alpha)
        {
            return (alpha.GetMean() - 1.0) * prob.Sum(Math.Log) - (prob.Count * alpha.GetMeanLogGamma(1.0) - alpha.GetMeanLogGamma(prob.Count)) ;
        }

        /// <summary>
        /// Evidence message for VMP
        /// </summary>
        /// <param name="prob">Constant value for 'prob'.</param>
        /// <param name="alpha">Incoming message for 'alpha'.</param>
        /// <returns>Average of the factor's log-value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>int p(alpha) log(factor(prob,alpha,K))</c>.
        /// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
        /// </para></remarks>
        [Quality(QualityBand.Experimental)]
        public static double AverageLogFactor([Proper] Dirichlet prob, ConjugateDirichlet alpha)
        {
            if (alpha.IsPointMass)
                return LogAverageFactor(prob, alpha.Point);
            return LogAverageFactor(prob, alpha.GetMode()); 
            //return (alpha.GetMean() - 1.0) * prob.GetMeanLog().Sum() - (prob.Dimension * alpha.GetMeanLogGamma(1.0) - alpha.GetMeanLogGamma(prob.Dimension));
        }


		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <param name="probMeanLog">Buffer for E[log(prob)]</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(prob) p(prob) log(factor(prob,alpha,K))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[Quality(QualityBand.Stable)]
		public static double AverageLogFactor([Proper] Dirichlet prob, double alpha, [Fresh] Vector probMeanLog)
		{
			double SumElogP = probMeanLog.Sum();
			double K = (double)probMeanLog.Count;
			return MMath.GammaLn(K*alpha) - K * MMath.GammaLn(alpha) + (alpha - 1.0) * SumElogP;
		}

		/// <summary>
		/// Returns the gradient and value of the KL divergence for this factor
		/// </summary>
		/// <param name="a2">Prior shape</param>
		/// <param name="b2">Prior rate</param>
		/// <param name="x">A vector of the variational posterior parameters. x[1]=log(shape), x[2]=log(rate)</param>
		/// <param name="SumElogP">Sum E[ log(prob_k) ]. Cached for efficiency</param>
		/// <param name="grad">Vector to fill with the gradient</param>
		/// <param name="K">Dimensionality</param>
		/// <returns>KL divergence</returns>
		private static double GradientAndValueAtPoint(double a2, double b2, Vector x, double SumElogP, Vector grad, double K)
		{
			double a = Math.Exp(x[0]);
			double b = Math.Exp(x[1]);
			double averageFactor =  GammaFromShapeAndRateOp.ELogGamma(Gamma.FromShapeAndRate(a, b / K));
			averageFactor -= K * GammaFromShapeAndRateOp.ELogGamma(Gamma.FromShapeAndRate(a, b));
			averageFactor += (a / b - 1.0) * SumElogP;
			double kl_value = Math.Log(b) - MMath.GammaLn(a) + (a - 1) * MMath.Digamma(a) - a // entropy
                - (a2 * Math.Log(b2) - MMath.GammaLn(a2) + (a2 - 1) * (MMath.Digamma(a) - Math.Log(b)) - b2 * a / b) // cross entropy
                - averageFactor; // factor
			if (double.IsInfinity(kl_value) || double.IsNaN(kl_value))
				throw new ApplicationException("KL divergence became ill-defined.");
			if (grad != null) {
				var gradS = Vector.Zero(2);
				gradS += GammaFromShapeAndRateOp.CalculateDerivatives(Gamma.FromShapeAndRate(a, b/K));
				gradS[1] = gradS[1] / K;
				gradS -= K * GammaFromShapeAndRateOp.CalculateDerivatives(Gamma.FromShapeAndRate(a, b));
				gradS[0] += SumElogP / b;
				gradS[1] -= SumElogP * a / (b * b);
				grad[0] = (a - 1.0) * MMath.Trigamma(a) - 1.0 // entropy
                    - (a2 - 1) * MMath.Trigamma(a) + b2 / b // cross term
                   - gradS[0]; // factor
				grad[0] *= a; // chain rule
				grad[1] = 1.0 / b // entropy
                    + (a2 - 1) / b - b2 * a / (b * b) // cross term
                   - gradS[1]; // factor
				grad[1] *= b; // chain rule
			}

			return kl_value;
		}

        /// <summary>
        /// Initialisation of buffer for E[log(prob)]
        /// </summary>
        /// <param name="prob">Incoming message from 'prob'</param>
        /// <returns>Correct size message.</returns>
        public static Vector ProbMeanLogInit([IgnoreDependency] Vector prob)
        {
            return Vector.Copy(prob);
        }

		/// <summary>
		/// Initialisation of buffer for E[log(prob)]
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'</param>
		/// <returns>Correct size message.</returns>
		public static Vector ProbMeanLogInit([IgnoreDependency] Dirichlet prob)
		{
			return Vector.Copy(prob.PseudoCount);
		}

		/// <summary>
		/// Buffer for E[log(prob)]
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'</param>
		/// <param name="result">Will be the returned value. </param>
		/// <returns>E[log(prob)]</returns>
		public static Vector ProbMeanLog(Dirichlet prob, Vector result)
		{
            prob.GetMeanLog(result);
            return result; 
		}

        /// <summary>
        /// Buffer for E[log(prob)]
        /// </summary>
        /// <param name="prob">Incoming message from 'prob'</param>
        /// <param name="result">Will be the returned value. </param>
        /// <returns>E[log(prob)]</returns>
        public static Vector ProbMeanLog(Vector prob, Vector result)
        {
            result.SetToFunction(prob, Math.Log); 
            return result;
        }

		/// <summary>
		/// VMP message to 'alpha'
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="alpha">Incoming message from 'alpha'.</param>
		/// <param name="probMeanLog">Buffer for E[log(prob)].</param>
		/// <param name="to_Alpha">Previous message sent to 'alpha'</param>
		/// <returns>Message to alpha</returns>
		/// <remarks><para>
		/// Optimal message calculated by minimising local KL divergence using LBFGS. 
		/// </para></remarks>
		[Quality(QualityBand.Experimental)]
		public static Gamma AlphaAverageLogarithm([Proper] Dirichlet prob, [SkipIfUniform] Gamma alpha, [Fresh] Vector probMeanLog, Gamma to_Alpha)
		{
			if (alpha.IsPointMass)
				return Gamma.Uniform();
			var s = new BFGS();
			int K = probMeanLog.Count;
			var prior = alpha / to_Alpha;
			int evalCounter = 0;
			s.MaximumStep = 20;
			s.MaximumIterations = 100;
			s.Epsilon = 1e-5;
			s.convergenceCriteria = BFGS.ConvergenceCriteria.Objective;
			double SumElogP = probMeanLog.Sum();
			var z = Vector.FromArray(new double[] { Math.Log(alpha.Shape), Math.Log(alpha.Rate) });
			double startingValue = GradientAndValueAtPoint(prior.Shape, prior.Rate, z, SumElogP, null, (double)K);
			FunctionEval f = delegate(Vector y, ref Vector grad) { evalCounter++; return GradientAndValueAtPoint(prior.Shape, prior.Rate, y, SumElogP, grad, (double)K); };
			//DerivativeChecker.CheckDerivatives(f, z); 
			z = s.Run(z, 1.0, f);
			var result = Gamma.FromShapeAndRate(Math.Exp(z[0]), Math.Exp(z[1]));
			result.SetToRatio(result, prior);
			double endValue = GradientAndValueAtPoint(prior.Shape, prior.Rate, z, SumElogP, null, (double)K);
			//Console.WriteLine("Went from {0} to {1} in {2} steps, {3} evals", startingValue, endValue, s.IterationsPerformed, evalCounter);
			if (startingValue < endValue)
				Console.WriteLine("Warning: BFGS resulted in an increased objective function");
			return result;
		}

        /// <summary>
        /// VMP message to 'alpha'
        /// </summary>
        /// <param name="prob">Constant value for prob 'prob'.</param>
        /// <param name="alpha">Incoming message from 'alpha'.</param>
        /// <param name="probMeanLog">Buffer for E[log(x)] values of 'prob'</param>
        /// <returns>Message to alpha</returns>
        /// <remarks><para>
        /// Optimal message calculated by minimising local KL divergence using LBFGS. 
        /// </para></remarks>
        [Quality(QualityBand.Experimental)]
        public static ConjugateDirichlet AlphaAverageLogarithm([Proper, SkipIfUniform] Dirichlet prob, [SkipIfUniform] ConjugateDirichlet alpha, [Fresh, SkipIfUniform] Vector probMeanLog)
        {
            // TODO: why is probMeanLog not set correctly? 
            if (alpha.IsPointMass)
                return ConjugateDirichlet.Uniform();
            var result = ConjugateDirichlet.FromNatural(0, -prob.GetMeanLog().Sum(), probMeanLog.Count, 1);
            //var result = ConjugateDirichlet.FromNatural(0, -probMeanLog.Sum(), probMeanLog.Count, 1);
            return result; 
        }

        /// <summary>
        /// VMP message to 'alpha'
        /// </summary>
        /// <param name="prob">Constant value for prob 'prob'.</param>
        /// <param name="alpha">Incoming message from 'alpha'.</param>
        /// <param name="probMeanLog">Buffer for E[log(x)] values of 'prob'</param>
        /// <returns>Message to alpha</returns>
        /// <remarks><para>
        /// Optimal message calculated by minimising local KL divergence using LBFGS. 
        /// </para></remarks>
        [Quality(QualityBand.Experimental)]
        public static ConjugateDirichlet AlphaAverageLogarithm([Proper] Vector prob, ConjugateDirichlet alpha, [Fresh, SkipIfUniform] Vector probMeanLog)
        {
            // TODO: why is probMeanLog not set correctly? 
            if (alpha.IsPointMass)
                return ConjugateDirichlet.Uniform();
            var result = ConjugateDirichlet.FromNatural(0, -prob.Sum(Math.Log), probMeanLog.Count, 1);
            return result;
        }

		/// <summary>
		/// VMP message to 'alpha'
		/// </summary>
		/// <param name="prob">Constant value for prob 'prob'.</param>
		/// <param name="alpha">Incoming message from 'alpha'.</param>
		/// <param name="to_Alpha">Previous message sent to 'alpha'</param>
		/// <returns>Message to alpha</returns>
		/// <remarks><para>
		/// Optimal message calculated by minimising local KL divergence using LBFGS. 
		/// </para></remarks>
		[Quality(QualityBand.Experimental)]
		public static Gamma AlphaAverageLogarithm(Vector prob, [SkipIfUniform] Gamma alpha, Gamma to_Alpha)
		{
			Dirichlet d = Dirichlet.PointMass(prob);
			return AlphaAverageLogarithm(d, alpha, d.GetMeanLog(), to_Alpha);
		}

		/// <summary>
		/// VMP message to 'prob'
		/// </summary>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is simply the definition of the factor since 
		/// alpha is fixed.
		/// </para></remarks>
		[Quality(QualityBand.Stable)]
		public static Dirichlet ProbAverageLogarithm([SkipIfUniform] double alpha, Dirichlet result)
		{
			result.PseudoCount.SetAllElementsTo(alpha);
            result.TotalCount = result.PseudoCount.Sum(); 
			return result;
		}

		[Skip]
		public static Dirichlet ProbAverageLogarithmInit(int K)
		{
			return Dirichlet.Uniform(K);
		}

		/// <summary>
		/// VMP message to 'prob'
		/// </summary>
		/// <param name="alpha">Incoming message from 'alpha'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'prob'.
		/// The formula is <c>exp(sum_(alpha) p(alpha) log(factor(prob,alpha,K)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="alpha"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Dirichlet ProbAverageLogarithm([SkipIfUniform] Gamma alpha, Dirichlet result)
		{
			result.PseudoCount.SetAllElementsTo(alpha.GetMean());
            result.TotalCount = result.PseudoCount.Sum(); 
			return result;
		}

        /// <summary>
        /// VMP message to 'prob'
        /// </summary>
        /// <param name="alpha">Incoming message from 'alpha'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
        /// <param name="result">Modified to contain the outgoing message</param>
        /// <returns><paramref name="result"/></returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'prob'.
        /// The formula is <c>exp(sum_(alpha) p(alpha) log(factor(prob,alpha,K)))</c>.
        /// </para></remarks>
        [Quality(QualityBand.Experimental)]
        public static Dirichlet ProbAverageLogarithm([SkipIfUniform] ConjugateDirichlet alpha, Dirichlet result)
        {
            if (!alpha.IsPointMass)
            {
                double mean, variance;
                alpha.SmartProposal(out mean, out variance);
                double alphaMode = Math.Exp(mean);
                if (double.IsNaN(alphaMode) || double.IsInfinity(alphaMode))
                    throw new ApplicationException("Nan message in ProbAverageLogarithm");
                //result.PseudoCount.SetAllElementsTo(alphaMode);
                result.PseudoCount.SetAllElementsTo(alpha.GetMean());
            }
            else
                result.PseudoCount.SetAllElementsTo(alpha.Point); 
            result.TotalCount = result.PseudoCount.Sum(); 
            return result;
        }


		// --------------------- EP not supported ------------------------------ 

		const string NotSupportedMessage = "Expectation Propagation does not currently support Dirichlet distributions with stochastic arguments";

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="alpha">Incoming message from 'alpha'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(prob,alpha) p(prob,alpha) factor(prob,alpha,K))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="alpha"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Dirichlet prob, [SkipIfUniform] Gamma alpha)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log factor(prob,alpha,K))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(IList<double> prob, double alpha)
		{
			int dim = prob.Count;
			double sum = MMath.GammaLn(dim*alpha) - dim*MMath.GammaLn(alpha);
			double alphaM1 = alpha-1;
			for (int i = 0; i < dim; i++) {
				sum += alphaM1*Math.Log(prob[i]);
			}
			return sum;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log factor(prob,alpha,K))</c>.
		/// </para></remarks>
		public static double LogEvidenceRatio(IList<double> prob, double alpha) { return LogAverageFactor(prob, alpha); }


		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log factor(prob,alpha,K))</c>.
		/// </para></remarks>
		public static double AverageLogFactor(IList<double> prob, double alpha) { return LogAverageFactor(prob, alpha); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log int_prob p(prob) factor(prob,alpha,K))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Dirichlet prob, double alpha)
		{
			int dim = prob.Dimension;
			double sum = MMath.GammaLn(dim*alpha) - dim*MMath.GammaLn(alpha) - prob.GetLogNormalizer() - MMath.GammaLn(prob.TotalCount+dim*alpha);
			for (int i = 0; i < dim; i++) {
				sum += MMath.GammaLn(prob.PseudoCount[i]+alpha);
			}
			return sum;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log int_prob p(prob) factor(prob,alpha,K))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Dirichlet prob, double alpha) { return 0.0; }

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="alpha">Constant value for 'alpha'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'prob' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet ProbAverageConditional([SkipIfUniform] double alpha, Dirichlet result)
		{
			result.PseudoCount.SetAllElementsTo(alpha);
            result.TotalCount = result.PseudoCount.Sum(); 
			return result;
		}

		[Skip]
		public static Dirichlet ProbAverageConditionalInit(int K)
		{
			return Dirichlet.Uniform(K);
		}

		/// <summary>
		/// EP message to 'alpha'
		/// </summary>
		/// <param name="alpha">Incoming message from 'alpha'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'alpha' as the random arguments are varied.
		/// The formula is <c>proj[p(alpha) sum_(prob) p(prob) factor(prob,alpha,K)]/p(alpha)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="alpha"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma AlphaAverageConditional([SkipIfUniform] Gamma alpha, [SkipIfUniform] Dirichlet prob, Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="alpha">Incoming message from 'alpha'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(alpha) p(alpha) factor(prob,alpha,K)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="alpha"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet ProbAverageConditional([SkipIfUniform] Gamma alpha, Dirichlet result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.DirichletFromMeanAndTotalCount"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "DirichletFromMeanAndTotalCount")]
	[Quality(QualityBand.Preview)]
	public static class DirichletOp
	{
		/// <summary>
		/// How much damping to use to prevent improper messages. Higher values result in more damping. 
		/// </summary>
		static public double damping = 0.0;

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Vector prob, Vector mean, double totalCount)
		{
			return (new Dirichlet(mean * totalCount)).GetLogProb(prob);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(prob,mean,totalCount) p(prob,mean,totalCount) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static double AverageLogFactor([SkipIfUniform] Dirichlet prob, [SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount)
		{
			double totalCountMean = totalCount.GetMean();
			Vector meanMean = mean.GetMean();
			Vector probMeanLog = prob.GetMeanLog();
			double sum = probMeanLog.Inner(meanMean, x => totalCountMean * x - 1.0);
			return sum + GammaFromShapeAndRateOp.ELogGamma(totalCount) - EvidenceMessageExpectations(mean, totalCount).Sum();
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(prob,totalCount) p(prob,totalCount) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static double AverageLogFactor([SkipIfUniform] Dirichlet prob, Vector mean, [SkipIfUniform] Gamma totalCount)
		{
			double totalCountMean = totalCount.GetMean();
			Vector probMeanLog = prob.GetMeanLog();
			double sum = GammaFromShapeAndRateOp.ELogGamma(totalCount);
			Gamma smk = new Gamma(totalCount);
			sum += probMeanLog.Inner(mean, x => totalCountMean * x - 1.0);
			sum += mean.Sum(x => { smk.Rate = totalCount.Rate / x; return -GammaFromShapeAndRateOp.ELogGamma(smk); });
			return sum;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(totalCount) p(totalCount) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Vector prob, Vector mean, [SkipIfUniform] Gamma totalCount)
		{
			return AverageLogFactor(Dirichlet.PointMass(prob), mean, totalCount);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(prob,mean) p(prob,mean) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor([SkipIfUniform] Dirichlet prob, [SkipIfUniform] Dirichlet mean, double totalCount)
		{
			return AverageLogFactor(prob, mean, Gamma.PointMass(totalCount));
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(prob) p(prob) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static double AverageLogFactor([SkipIfUniform] Dirichlet prob, Vector mean, double totalCount)
		{
			return AverageLogFactor(prob, mean, Gamma.PointMass(totalCount));
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean,totalCount) p(mean,totalCount) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Vector prob, [SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount)
		{
			return AverageLogFactor(Dirichlet.PointMass(prob), mean, totalCount);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean) p(mean) log(factor(prob,mean,totalCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Vector prob, [SkipIfUniform] Dirichlet mean, double totalCount)
		{
			return AverageLogFactor(Dirichlet.PointMass(prob), mean, Gamma.PointMass(totalCount));
		}

		/// <summary>
		/// VMP message to 'prob'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'prob'.
		/// The formula is <c>exp(sum_(mean,totalCount) p(mean,totalCount) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Dirichlet ProbAverageLogarithm([SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount)
		{
			return new Dirichlet(totalCount.GetMean() * mean.GetMean());
		}

		/// <summary>
		/// VMP message to 'prob'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'prob'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static Dirichlet ProbAverageLogarithm([SkipIfUniform] Dirichlet mean, double totalCount)
		{
			return new Dirichlet(totalCount * mean.GetMean());
		}

		/// <summary>
		/// VMP message to 'prob'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'prob'.
		/// The formula is <c>exp(sum_(totalCount) p(totalCount) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Dirichlet ProbAverageLogarithm(Vector mean, [SkipIfUniform] Gamma totalCount)
		{
			return new Dirichlet(totalCount.GetMean() * mean);
		}

		/// <summary>
		/// VMP message to 'prob'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'prob' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet ProbAverageLogarithm(Vector mean, double totalCount)
		{
			return new Dirichlet(totalCount * mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_mean">Previous outgoing message to 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(totalCount) p(totalCount) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Dirichlet MeanAverageLogarithm([Proper] Dirichlet mean, [Proper] Gamma totalCount, Vector prob, Dirichlet to_mean)
		{
			return MeanAverageLogarithm(mean, totalCount, Dirichlet.PointMass(prob), to_mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_mean">Previous outgoing message to 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static Dirichlet MeanAverageLogarithm([Proper] Dirichlet mean, double totalCount, Vector prob, Dirichlet to_mean)
		{
			return MeanAverageLogarithm(mean, Gamma.PointMass(totalCount), Dirichlet.PointMass(prob), to_mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_mean">Previous outgoing message to 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(prob) p(prob) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Dirichlet MeanAverageLogarithm([Proper] Dirichlet mean, double totalCount, [SkipIfUniform] Dirichlet prob, Dirichlet to_mean)
		{
			return MeanAverageLogarithm(mean, Gamma.PointMass(totalCount), prob, to_mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_mean">Previous outgoing message to 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(totalCount,prob) p(totalCount,prob) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Dirichlet MeanAverageLogarithm([Proper] Dirichlet mean, [Proper] Gamma totalCount, [SkipIfUniform] Dirichlet prob, Dirichlet to_mean)
		{
			Vector gradS = CalculateGradientForMean(mean.PseudoCount, totalCount, prob.GetMeanLog());
			// Project onto Dirichlet, efficient matrix inversion (see TM's Dirichlet fitting paper)
			int K = mean.Dimension;
			Vector q = Vector.Zero(K);
			double gOverQ = 0, OneOverQ = 0;
			for (int k = 0; k < K; k++) {
				q[k] = MMath.Trigamma(mean.PseudoCount[k]);
				gOverQ += gradS[k] / q[k];
				OneOverQ += 1 / q[k];
			}
			double z = -MMath.Trigamma(mean.TotalCount);
			double b = gOverQ / (1 / z + OneOverQ);
			// Create new approximation and damp

			if (damping == 0.0) {
				to_mean.PseudoCount.SetToFunction(gradS, q, (x, y) => ((x - b) / y) + 1.0);
				return to_mean;
			} else {
				var old_msg = (Dirichlet)to_mean.Clone();
				to_mean.PseudoCount.SetToFunction(gradS, q, (x, y) => ((x - b) / y) + 1.0);
				return (to_mean ^ (1 - damping)) * (old_msg ^ damping);
			}
		}

		/// <summary>
		/// Helper function to calculate gradient of the KL divergence with respect to the mean of the Dirichlet. 
		/// </summary>
		/// <param name="meanPseudoCount">Pseudocount vector of the incoming message from 'mean'</param>
		/// <param name="totalCount">Incoming message from 'totalCount'</param>
		/// <param name="meanLogProb">E[log(prob)]</param>
		/// <returns>Gradient of the KL divergence with respect to the mean of the Dirichlet</returns>
		internal static Vector CalculateGradientForMean(Vector meanPseudoCount, Gamma totalCount, Vector meanLogProb)
		{
			// Compute required integrals
			double[] EELogGamma;
			double[] EELogMLogGamma;
			double[] EELogOneMinusMLogGamma;
			MeanMessageExpectations(
					meanPseudoCount,
					totalCount,
					out EELogGamma,
					out EELogMLogGamma,
					out EELogOneMinusMLogGamma);

			// Calculate gradients of ELogGamma(sm)
			int K = meanPseudoCount.Count;
			double meanTotalCount = meanPseudoCount.Sum();
			Vector ELogM = Vector.Zero(K);
			Vector B = Vector.Zero(K);
			Vector A = Vector.Zero(K);
			ELogM.SetToFunction(meanPseudoCount, MMath.Digamma);
			ELogM.SetToDifference(ELogM, MMath.Digamma(meanTotalCount));
			for (int k = 0; k < K; k++) {
				A[k] = EELogMLogGamma[k] - ELogM[k] * EELogGamma[k];
				double ELogOneMinusM = MMath.Digamma(meanTotalCount - meanPseudoCount[k])
                    - MMath.Digamma(meanTotalCount);
				B[k] = EELogOneMinusMLogGamma[k] - ELogOneMinusM * EELogGamma[k];
			}
			Vector gradC = A - B + B.Sum();
			// Calculate gradients of analytic part
			double sum = 0;
			for (int k = 0; k < K; k++)
				sum += meanPseudoCount[k] * meanLogProb[k];
			Vector gradS = Vector.Constant(K, -sum / (meanTotalCount * meanTotalCount));
			for (int k = 0; k < K; k++)
				gradS[k] += meanLogProb[k] / meanTotalCount;
			gradS *= totalCount.GetMean();
			gradS -= gradC;
			return gradS;
		}

		/// <summary>
		/// VMP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_totalCount">Previous outgoing message to 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'totalCount' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'totalCount' conditioned on the given values.
		/// </para>
		/// <para>
		/// The outgoing message here would not be Dirichlet distributed, so we use Nonconjugate VMP, which
		/// sends the approximate factor ensuring the gradient of the KL wrt to the variational parameters match. 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Gamma TotalCountAverageLogarithm(Vector mean, [Proper] Gamma totalCount, Vector prob, Gamma to_totalCount)
		{
			return TotalCountAverageLogarithm(mean, totalCount, Dirichlet.PointMass(prob), to_totalCount);
		}

		/// <summary>
		/// VMP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_totalCount">Previous outgoing message to 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'totalCount' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'totalCount'.
		/// The formula is <c>exp(sum_(prob) p(prob) log(factor(prob,mean,totalCount)))</c>.
		/// </para>
		/// <para>
		/// The outgoing message here would not be Dirichlet distributed, so we use Nonconjugate VMP, which
		/// sends the approximate factor ensuring the gradient of the KL wrt to the variational parameters match. 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Gamma TotalCountAverageLogarithm(Vector mean, [Proper] Gamma totalCount, [SkipIfUniform] Dirichlet prob, Gamma to_totalCount)
		{
			double at = totalCount.Shape;
			double bt = totalCount.Rate;
			// Find required expectations using quadrature
			Vector gradElogGamma = GammaFromShapeAndRateOp.CalculateDerivatives(totalCount);
			Vector gradS = gradElogGamma;
			Gamma smk = new Gamma(totalCount);
			for (int k = 0; k < mean.Count; k++) {
				smk.Rate = totalCount.Rate / mean[k];
				gradS -= GammaFromShapeAndRateOp.CalculateDerivatives(smk);
			}
			// Analytic 
			double c = mean.Inner(prob.GetMeanLog());
			gradS[0] += c / bt;
			gradS[1] -= c * at / (bt * bt);
			Matrix mat = new Matrix(2, 2);
			mat[0, 0] = MMath.Trigamma(at);
			mat[1, 0] = mat[0, 1] = -1 / bt;
			mat[1, 1] = at / (bt * bt);
			Vector v = GammaFromShapeAndRateOp.twoByTwoInverse(mat) * gradS;
			Gamma approximateFactor = Gamma.FromShapeAndRate(v[0] + 1, v[1]);

			if (damping == 0.0)
				return approximateFactor;
			else
				return (approximateFactor ^ (1 - damping)) * (to_totalCount ^ damping);
		}

		/// <summary>
		/// VMP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_totalCount">Previous outgoing message to 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'totalCount' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'totalCount'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(prob,mean,totalCount)))</c>.
		/// </para><para>
		/// The outgoing message here would not be Dirichlet distributed, so we use Nonconjugate VMP, which
		/// sends the approximate factor ensuring the gradient of the KL wrt to the variational parameters match. 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Gamma TotalCountAverageLogarithm([Proper] Dirichlet mean, [Proper] Gamma totalCount, Vector prob, Gamma to_totalCount)
		{
			return TotalCountAverageLogarithm(mean, totalCount, Dirichlet.PointMass(prob), to_totalCount);
		}

		/// <summary>
		/// VMP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_totalCount">Previous outgoing message to 'totalCount'.</param>
		/// <remarks><para>
		/// The outgoing message here would not be Dirichlet distributed, so we use Nonconjugate VMP, which
		/// sends the approximate factor ensuring the gradient of the KL wrt to the variational parameters match. 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Gamma TotalCountAverageLogarithm([Proper] Dirichlet mean, [Proper] Gamma totalCount, [SkipIfUniform] Dirichlet prob, Gamma to_totalCount)
		{
			Gamma approximateFactor = TotalCountAverageLogarithmHelper(
					mean.PseudoCount,
					totalCount,
					prob.GetMeanLog());
			double damping = 0.9;
			if (damping == 0.0)
				return approximateFactor;
			else
				return (approximateFactor ^ (1 - damping)) * (to_totalCount ^ damping);
		}

		/// <summary>
		/// VMP message to 'totalCount'. This functionality is separated out to allow use by BetaOp. 
		/// </summary>
		/// <param name="meanPseudoCount">Pseudocount of incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="meanLogProb">E[log(prob)] from incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <remarks><para>
		/// The outgoing message here would not be Dirichlet distributed, so we use Nonconjugate VMP, which
		/// sends the approximate factor ensuring the gradient of the KL wrt to the variational parameters match. 
		/// </para></remarks>
		internal static Gamma TotalCountAverageLogarithmHelper(Vector meanPseudoCount, Gamma totalCount, Vector meanLogProb)
		{
			double[] EELogGamma;
			double[] EELogSLogGamma;
			double[] EEMSDigamma;
			// 2D quadrature
			TotalCountMessageExpectations(
					meanPseudoCount,
					totalCount,
					out EELogGamma,
					out EELogSLogGamma,
					out EEMSDigamma);
			double at = totalCount.Shape;
			double bt = totalCount.Rate;
			// Find required expectations using quadrature
			Vector gradElogGamma = GammaFromShapeAndRateOp.CalculateDerivatives(totalCount);
			Vector gradS = gradElogGamma;
			Vector EM = Vector.Zero(meanPseudoCount.Count);
			EM.SetToProduct(meanPseudoCount, 1.0 / meanPseudoCount.Sum());
			double c = 0;
			for (int k = 0; k < meanPseudoCount.Count; k++) {
				gradS[0] -= EELogSLogGamma[k] - totalCount.GetMeanLog() * EELogGamma[k];
				gradS[1] -= -EEMSDigamma[k] / bt;
				c += EM[k] * meanLogProb[k];
			}
			// Analytic 
			gradS[0] += c / bt;
			gradS[1] -= c * at / (bt * bt);
			Matrix mat = new Matrix(2, 2);
			mat[0, 0] = MMath.Trigamma(at);
			mat[1, 0] = mat[0, 1] = -1 / bt;
			mat[1, 1] = at / (bt * bt);
			Vector v = GammaFromShapeAndRateOp.twoByTwoInverse(mat) * gradS;
			return Gamma.FromShapeAndRate(v[0] + 1, v[1]);
		}

		/// <summary>
		/// Perform the quadrature required for the Nonconjugate VMP message to 'mean'
		/// </summary>
		/// <param name="meanQPseudoCount">Incoming message from 'mean'.</param>
		/// <param name="totalCountQ">Incoming message from 'totalCount'.</param>
		/// <param name="EELogGamma">Array to be filled with E[LogGamma(s*m_k)].</param>
		/// <param name="EELogMLogGamma">Array to be filled with E[Log(m_k)*LogGamma(s*m_k)].</param>
		/// <param name="EELogOneMinusMLogGamma">Array to be filled with E[Log(1-m_k)*LogGamma(s*m_k)].</param>
		/// <remarks><para>
		/// All three arrays are calculated simultaneously for efficiency. The quadrature over 
		/// 'totalCount' (which is Gamma-distributed) is peformed by a change of variable x=log(s)
		/// followed by Gauss-Hermite quadrature. The quadrature over m is performed using 
		/// Gauss-Legendre. 
		/// </para></remarks>
		public static void MeanMessageExpectations(
				Vector meanQPseudoCount,
				Gamma totalCountQ,
				out double[] EELogGamma,
				out double[] EELogMLogGamma,
				out double[] EELogOneMinusMLogGamma)
		{
			// Get shape and scale of the distribution
			double at, bt;
			at = totalCountQ.Shape;
			bt = totalCountQ.Rate;

			// Mean in the transformed domain
			double ELogS = totalCountQ.GetMeanLog();
			// Laplace approximation of variance in transformed domain 
			double proposalVariance = 1 / at;

			// Quadrature coefficient
			int nt = 32;
			Vector nodes = Vector.Zero(nt);
			Vector weights = Vector.Zero(nt);
			Vector expx = Vector.Zero(nt);
			if (!totalCountQ.IsPointMass) {
				Quadrature.GaussianNodesAndWeights(ELogS, proposalVariance, nodes, weights);
				// Precompute weights for each m slice
				for (int i = 0; i < nt; i++) {
					double x = nodes[i];
					expx[i] = Math.Exp(x);
					double p = at * x - bt * expx[i] - Gaussian.GetLogProb(x, ELogS, proposalVariance);
					weights[i] *= Math.Exp(p);
				}
			}

			int nm = 20;
			Vector mnodes = Vector.Zero(nm);
			Vector mweight = Vector.Zero(nm);
			Quadrature.UniformNodesAndWeights(0, 1, mnodes, mweight);
			int K = meanQPseudoCount.Count;
			Vector[] mweights = new Vector[K];
			Beta[] mkDist = new Beta[K];
			EELogGamma = new double[K];
			EELogMLogGamma = new double[K];
			EELogOneMinusMLogGamma = new double[K];
			double meanQTotalCount = meanQPseudoCount.Sum();
			for (int i = 0; i < K; i++) {
				mweights[i] = Vector.Copy(mweight);
				mkDist[i] = new Beta(meanQPseudoCount[i], meanQTotalCount - meanQPseudoCount[i]);
				EELogGamma[i] = 0;
				EELogMLogGamma[i] = 0;
				EELogOneMinusMLogGamma[i] = 0;
			}

			double ES = totalCountQ.GetMean();
			double ESLogS = ELogS * ES + 1 / bt;

			for (int j = 0; j < nm; j++) {
				double m = mnodes[j];
				double ELogGamma = 0;
				if (totalCountQ.IsPointMass)
					ELogGamma = MMath.GammaLn(m * totalCountQ.Point);
				else {
					// Calculate expectations in x=log(s) space using Gauss-Hermite quadrature
					for (int i = 0; i < nt; i++)
						ELogGamma += weights[i] * (MMath.GammaLn(m * expx[i]) + nodes[i]);
					// Normalise and add removed components
					double normalisation = Math.Pow(bt, at) / MMath.Gamma(at);
					ELogGamma = normalisation * ELogGamma - ELogS;
				}

				double EELogMLogGammaTemp = Math.Log(m) * (ELogGamma + ELogS + Math.Log(m));
				double EELogOneMinusMLogGammaTemp = Math.Log(1 - m) *
                    (ELogGamma - (.5 * Math.Log(2 * Math.PI) - .5 * ELogS
                    - .5 * Math.Log(m) + m * ESLogS + ES * m * Math.Log(m) - ES * m));
				for (int i = 0; i < K; i++) {
					mweights[i][j] *= Math.Exp(mkDist[i].GetLogProb(m));
					EELogGamma[i] += mweights[i][j] * ELogGamma;
					EELogMLogGamma[i] += mweights[i][j] * EELogMLogGammaTemp;
					EELogOneMinusMLogGamma[i] += mweights[i][j] * EELogOneMinusMLogGammaTemp;
				}
			}
			for (int i = 0; i < K; i++)
				AddAnalyticComponent(
						mkDist[i],
						ELogS,
						ES,
						ESLogS,
						ref EELogMLogGamma[i],
						ref EELogOneMinusMLogGamma[i]);
		}

		// Helper function to add the removed parts back (see note)
		private static void AddAnalyticComponent(
				Beta meanQ,
				double ELogS,
				double ES,
				double ESLogS,
				ref double EELogMLogGamma,
				ref double EELogOneMinusMLogGamma)
		{
			double ELogM, ELogOneMinusM;
			meanQ.GetMeanLogs(out ELogM, out ELogOneMinusM);
			double ELogMSquared = ELogM * ELogM
                + MMath.Trigamma(meanQ.TrueCount) - MMath.Trigamma(meanQ.TotalCount);
			EELogMLogGamma -= ELogS * ELogM + ELogMSquared;
			double Em = meanQ.GetMean();
			double am = meanQ.TrueCount;
			double bm = meanQ.FalseCount;
			double EmlogOneMinusM = Em * (MMath.Digamma(bm) - MMath.Digamma(am + bm + 1));
			double EmlogmlogOneMinusM = Em * ((MMath.Digamma(bm) - MMath.Digamma(am + bm + 1)) * (MMath.Digamma(am + 1)
                - MMath.Digamma(am + bm + 1)) - MMath.Trigamma(am + bm + 1));
			double ELogMLogOneMinusM = ELogM * ELogOneMinusM - MMath.Trigamma(meanQ.TotalCount);
			EELogOneMinusMLogGamma += .5 * Math.Log(2 * Math.PI) * ELogOneMinusM - .5 * ELogS * ELogOneMinusM
                - .5 * ELogMLogOneMinusM + EmlogOneMinusM * ESLogS
                + ES * EmlogmlogOneMinusM - ES * EmlogOneMinusM;
		}

		/// <summary>
		/// Perform the quadrature required for the Nonconjugate VMP message to 'totalCount'
		/// </summary>
		/// <param name="meanQPseudoCount">Incoming message from 'mean'.</param>
		/// <param name="totalCountQ">Incoming message from 'totalCount'.</param>
		/// <param name="EELogGamma">Array to be filled with E[LogGamma(s*m_k)].</param>
		/// <param name="EELogSLogGamma">Array to be filled with E[Log(s)*LogGamma(s*m_k)].</param>
		/// <param name="EEMSDigamma">Array to be filled with E[s*m_k*Digamma(s*m_k)].</param>
		/// <remarks><para>
		/// All three arrays are calculated simultaneously for efficiency. The quadrature over 
		/// 'totalCount' (which is Gamma-distributed) is peformed by a change of variable x=log(s)
		/// followed by Gauss-Hermite quadrature. The quadrature over m is performed using 
		/// Gauss-Legendre. 
		/// </para></remarks>
		public static void TotalCountMessageExpectations(
				Vector meanQPseudoCount,
				Gamma totalCountQ,
				out double[] EELogGamma,
				out double[] EELogSLogGamma,
				out double[] EEMSDigamma)
		{
			// Get shape and rate of the distribution
			double at = totalCountQ.Shape, bt = totalCountQ.Rate;

			// Mean in the transformed domain
			double proposalMean = totalCountQ.GetMeanLog();
			// Laplace approximation of variance in transformed domain 
			double proposalVariance = 1 / at;

			// Quadrature coefficient
			int nt = 32;
			Vector nodes = Vector.Zero(nt);
			Vector weights = Vector.Zero(nt);
			Vector expx = Vector.Zero(nt);
			if (!totalCountQ.IsPointMass) {
				Quadrature.GaussianNodesAndWeights(proposalMean, proposalVariance, nodes, weights);
				// Precompute weights for each m slice
				for (int i = 0; i < nt; i++) {
					double x = nodes[i];
					expx[i] = Math.Exp(x);
					double p = at * x - bt * expx[i] - Gaussian.GetLogProb(x, proposalMean, proposalVariance);
					weights[i] *= Math.Exp(p);
				}
			}

			int nm = 20;
			Vector mnodes = Vector.Zero(nm);
			Vector mweight = Vector.Zero(nm);
			Quadrature.UniformNodesAndWeights(0, 1, mnodes, mweight);
			int K = meanQPseudoCount.Count;
			Vector[] mweights = new Vector[K];
			Beta[] mkDist = new Beta[K];
			EELogGamma = new double[K];
			EELogSLogGamma = new double[K];
			EEMSDigamma = new double[K];
			double meanQTotalCount = meanQPseudoCount.Sum();
			for (int i = 0; i < K; i++) {
				mweights[i] = Vector.Copy(mweight);
				mkDist[i] = new Beta(meanQPseudoCount[i], meanQTotalCount - meanQPseudoCount[i]);
				EELogGamma[i] = 0;
				EELogSLogGamma[i] = 0;
				EEMSDigamma[i] = 0;
			}

			for (int j = 0; j < nm; j++) {
				double m = mnodes[j];
				double ESDigamma = 0;
				double ELogGamma = 0;
				double ELogSLogGamma = 0;
				if (totalCountQ.IsPointMass) {
					ESDigamma = totalCountQ.Point * MMath.Digamma(m * totalCountQ.Point);
					ELogGamma = MMath.GammaLn(m * totalCountQ.Point);
					ELogSLogGamma = Math.Log(totalCountQ.Point) * ELogGamma;
				} else {
					// Calculate expectations in x=log(s) space using Gauss-Hermite quadrature
					for (int i = 0; i < nt; i++) {
						double x = nodes[i];
						ELogGamma += weights[i] * (MMath.GammaLn(m * expx[i]) + x);
						ESDigamma += weights[i] * (expx[i] * MMath.Digamma(m * expx[i]) + 1);
						ELogSLogGamma += weights[i] * (x * MMath.GammaLn(m * expx[i]) + x * x + x * Math.Log(m));
					}
					// Normalise and add removed components
					double normalisation = Math.Pow(bt, at) / MMath.Gamma(at);
					ELogGamma = normalisation * ELogGamma - proposalMean;
					ELogSLogGamma = normalisation * ELogSLogGamma
                    - (MMath.Trigamma(at) + proposalMean * proposalMean + Math.Log(m) * proposalMean);
					ESDigamma = normalisation * ESDigamma - 1;
				}
				for (int i = 0; i < K; i++) {
					mweights[i][j] *= Math.Exp(mkDist[i].GetLogProb(m));
					EELogGamma[i] += mweights[i][j] * ELogGamma;
					EELogSLogGamma[i] += mweights[i][j] * ELogSLogGamma;
					EEMSDigamma[i] += mweights[i][j] * m * ESDigamma;
				}
			}
		}

		/// <summary>
		/// Perform the quadrature required for the VMP evidence message
		/// </summary>
		/// <param name="meanQ">Incoming message from m='mean'.</param>
		/// <param name="totalCountQ">Incoming message from s='totalCount'.</param>
		/// <returns>Vector of E[ LogGamma(s*m_k)].</returns>
		/// <remarks><para>
		/// The quadrature over 'totalCount' (which is Gamma-distributed) is 
		/// peformed by a change of variable x=log(s) followed by Gauss-Hermite 
		/// quadrature. The quadrature over m is performed using Gauss-Legendre. 
		/// </para></remarks>
		public static Vector EvidenceMessageExpectations(
				Dirichlet meanQ,
				Gamma totalCountQ)
		{
			// Get shape and scale of the distribution
			double at, bt;
			totalCountQ.GetShapeAndScale(out at, out bt);
			bt = 1 / bt; // want rate not scale

			// Mean in the transformed domain
			double proposalMean = totalCountQ.GetMeanLog();
			// Laplace approximation of variance in transformed domain 
			double proposalVariance = 1 / at;

			// Quadrature coefficient
			int nt = 32;
			Vector nodes = Vector.Zero(nt);
			Vector weights = Vector.Zero(nt);
			Vector expx = Vector.Zero(nt);
			if (!totalCountQ.IsPointMass) {
				Quadrature.GaussianNodesAndWeights(proposalMean, proposalVariance, nodes, weights);
				// Precompute weights for each m slice
				for (int i = 0; i < nt; i++) {
					double x = nodes[i];
					expx[i] = Math.Exp(x);
					double p = at * x - bt * expx[i] - Gaussian.GetLogProb(x, proposalMean, proposalVariance);
					weights[i] *= Math.Exp(p);
				}
			}

			int nm = 20;
			Vector mnodes = Vector.Zero(nm);
			Vector mweight = Vector.Zero(nm);
			Quadrature.UniformNodesAndWeights(0, 1, mnodes, mweight);
			int K = meanQ.Dimension;
			Vector[] mweights = new Vector[K];
			Beta[] mkDist = new Beta[K];
			double[] EELogGamma = new double[K];
			for (int i = 0; i < K; i++) {
				mweights[i] = Vector.Copy(mweight);
				mkDist[i] = new Beta(meanQ.PseudoCount[i], meanQ.TotalCount - meanQ.PseudoCount[i]);
				EELogGamma[i] = 0;
			}

			for (int j = 0; j < nm; j++) {
				double m = mnodes[j];
				double ELogGamma = 0;
				if (totalCountQ.IsPointMass)
					ELogGamma = MMath.GammaLn(m * totalCountQ.Point);
				else {
					// Calculate expectations in x=log(s) space using Gauss-Hermite quadrature
					for (int i = 0; i < nt; i++) {
						double x = nodes[i];
						ELogGamma += weights[i] * (MMath.GammaLn(m * expx[i]) + x);
					}
					// Normalise and add removed components
					double normalisation = Math.Pow(bt, at) / MMath.Gamma(at);
					ELogGamma = normalisation * ELogGamma - proposalMean;
				}
				for (int i = 0; i < K; i++) {
					mweights[i][j] *= Math.Exp(mkDist[i].GetLogProb(m));
					EELogGamma[i] += mweights[i][j] * ELogGamma;
				}
			}
			return Vector.FromArray(EELogGamma);
		}


		// ------------------------------------ EP ----------------------------------

		const string NotSupportedMessage = "Expectation Propagation does not currently support Dirichlet distributions with stochastic arguments";

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector prob, Vector mean, double totalCount)
		{
			var temp = mean.Clone();
			mean.SetToProduct(mean, totalCount);
			var d = new Dirichlet(temp);
			return d.GetLogProb(prob);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector prob, Vector mean, double totalCount) { return LogAverageFactor(prob, mean, totalCount); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="to_prob">Previous outgoing message to 'prob'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(prob) p(prob) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Dirichlet prob, Vector mean, Dirichlet to_prob, double totalCount)
		{
			return to_prob.GetLogAverageOf(prob);
		}

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'prob' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet ProbAverageConditional(Dirichlet prob, Vector mean, double totalCount, Dirichlet result)
		{
			result.PseudoCount.SetToProduct(mean, totalCount);
			return result;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(prob,mean,totalCount) p(prob,mean,totalCount) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Dirichlet prob, [SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(prob,totalCount) p(prob,totalCount) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Dirichlet prob, Vector mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(totalCount) p(totalCount) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor(Vector prob, Vector mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(prob,mean) p(prob,mean) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Dirichlet prob, [SkipIfUniform] Dirichlet mean, double totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean,totalCount) p(mean,totalCount) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor(Vector prob, [SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean) p(mean) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor(Vector prob, [SkipIfUniform] Dirichlet mean, double totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(mean,totalCount) p(mean,totalCount) factor(prob,mean,totalCount)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet ProbAverageConditional([SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>The outgoing EP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(mean) p(mean) factor(prob,mean,totalCount)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet ProbAverageConditional([SkipIfUniform] Dirichlet mean, double totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(totalCount) p(totalCount) factor(prob,mean,totalCount)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet ProbAverageConditional(Vector mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>The outgoing EP message to the 'prob' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'prob' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet ProbAverageConditional(Vector mean, double totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(totalCount) p(totalCount) factor(prob,mean,totalCount)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet MeanAverageConditional([SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount, Vector prob, [SkipIfUniform] Dirichlet result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet MeanAverageConditional([SkipIfUniform] Dirichlet mean, double totalCount, Vector prob, [SkipIfUniform] Dirichlet result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(prob) p(prob) factor(prob,mean,totalCount)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet MeanAverageConditional([SkipIfUniform] Dirichlet mean, double totalCount, [SkipIfUniform] Dirichlet prob, [SkipIfUniform] Dirichlet result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(totalCount,prob) p(totalCount,prob) factor(prob,mean,totalCount)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Dirichlet MeanAverageConditional([SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount, [SkipIfUniform] Dirichlet prob, [SkipIfUniform] Dirichlet result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'totalCount' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma TotalCountAverageConditional(Vector mean, [SkipIfUniform] Gamma totalCount, Vector prob, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'totalCount' as the random arguments are varied.
		/// The formula is <c>proj[p(totalCount) sum_(prob) p(prob) factor(prob,mean,totalCount)]/p(totalCount)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma TotalCountAverageConditional(Vector mean, [SkipIfUniform] Gamma totalCount, [SkipIfUniform] Dirichlet prob, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'totalCount' as the random arguments are varied.
		/// The formula is <c>proj[p(totalCount) sum_(mean) p(mean) factor(prob,mean,totalCount)]/p(totalCount)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma TotalCountAverageConditional([SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount, Vector prob, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'totalCount' as the random arguments are varied.
		/// The formula is <c>proj[p(totalCount) sum_(mean,prob) p(mean,prob) factor(prob,mean,totalCount)]/p(totalCount)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma TotalCountAverageConditional([SkipIfUniform] Dirichlet mean, [SkipIfUniform] Gamma totalCount, [SkipIfUniform] Dirichlet prob, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
	}
}
