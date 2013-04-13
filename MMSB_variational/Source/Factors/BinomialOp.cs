// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Binomial"/>, given random arguments to the function.
	/// </summary>
	/// <remarks>The factor is f(sample,n,p) = choose(n,sample) p^sample (1-p)^(n-sample)</remarks>
	[FactorMethod(typeof(Factor), "Binomial", typeof(int), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class BinomialOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete sample, [Fresh] Discrete to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, double p, int trialCount)
		{
			return MMath.ChooseLn(trialCount, sample) + sample*Math.Log(p) + (trialCount-sample)*Math.Log(1-p);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(trialCount) p(trialCount) factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, double p, Discrete trialCount)
		{
			double logZ = Double.NegativeInfinity;
			for (int n = 0; n < trialCount.Dimension; n++)
			{
				logZ = MMath.LogSumExp(logZ, trialCount.GetLogProb(n) + LogAverageFactor(sample, p, n));
			}
			return logZ;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(p,trialCount) p(p,trialCount) factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, Beta p, Discrete trialCount)
		{
			double logZ = Double.NegativeInfinity;
			for (int n = 0; n < trialCount.Dimension; n++)
			{
				logZ = MMath.LogSumExp(logZ, trialCount.GetLogProb(n) + LogAverageFactor(sample, p, n));
			}
			return logZ;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(p) p(p) factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, Beta p, int trialCount)
		{
			double a = p.TrueCount;
			double b = p.FalseCount;
			return MMath.ChooseLn(trialCount, sample) +
				MMath.GammaLn(a+sample)-MMath.GammaLn(a) +
				MMath.GammaLn(b+trialCount-sample)-MMath.GammaLn(b) +
				MMath.GammaLn(a+b)-MMath.GammaLn(a+b+trialCount);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,trialCount,p) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Discrete sample) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(trialCount) p(trialCount) factor(sample,trialCount,p))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sample, double p, Discrete trialCount)
		{
			return LogAverageFactor(sample, p, trialCount);
		}
		public static double LogEvidenceRatio(int sample, double p, int trialCount)
		{
			return LogAverageFactor(sample, p, trialCount);
		}
		public static double LogEvidenceRatio(int sample, Beta p, int trialCount)
		{
			return LogAverageFactor(sample, p, trialCount);
		}

		[Skip]
		public static Discrete SampleAverageConditionalInit(int trialCount)
		{
			return Discrete.Uniform(trialCount+1);
		}
		[Skip]
		public static Discrete SampleAverageConditionalInit([IgnoreDependency] Discrete trialCount)
		{
			return Discrete.Uniform(trialCount.Dimension);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(double p, int trialCount, Discrete result)
		{
			// result must range from 0 to n
			if (result.Dimension < trialCount+1) throw new ArgumentException("result.Dimension ("+result.Dimension+") < n+1 ("+trialCount+"+1)");
			Vector probs = result.GetWorkspace();
			double logp = Math.Log(p);
			double log1minusp = Math.Log(1-p);
			probs.SetAllElementsTo(0.0);
			for (int k = 0; k <= trialCount; k++)
				probs[k] = Math.Exp(MMath.ChooseLn(trialCount, k) + k*logp + (trialCount-k)*log1minusp);
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(p) p(p) factor(sample,trialCount,p)]/p(sample)</c>.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(Beta p, int trialCount, Discrete result)
		{
			if (p.IsPointMass) return SampleAverageConditional(p.Point, trialCount, result);
			// result must range from 0 to n
			if (result.Dimension < trialCount+1) throw new ArgumentException("result.Dimension ("+result.Dimension+") < n+1 ("+trialCount+"+1)");
			Vector probs = result.GetWorkspace();
			double a = p.TrueCount;
			double b = p.FalseCount;
			probs.SetAllElementsTo(0.0);
			double max = double.NegativeInfinity;
			for (int k = 0; k <= trialCount; k++) {
				double logProb = MMath.GammaLn(a+k)-MMath.GammaLn(1+k) + MMath.GammaLn(b+trialCount-k)-MMath.GammaLn(1+trialCount-k);
				probs[k] = logProb;
				if (logProb > max) max = logProb;
			}
			for (int k = 0; k <= trialCount; k++) {
				probs[k] = Math.Exp(probs[k] - max);
			}
			result.SetProbs(probs);
			return result;
		}

		public static Poisson SampleAverageConditional(double p, Poisson trialCount)
		{
			if(trialCount.Precision != 1) throw new NotImplementedException("precision != 1 not implemented");
			// p(sample=k) = sum(n>=k) lambda^n/n! exp(-lambda) nchoosek(n,k) p^k (1-p)^(n-k)
			//             = (p/(1-p))^k 1/k! exp(-lambda) sum_(n>=k) 1/(n-k)! (1-p)^n lambda^n
			//             = (p/(1-p))^k 1/k! exp(-lambda) (1-p)^k lambda^k exp((1-p) lambda)
			//             = p^k lambda^k 1/k! exp(-p lambda)
			return new Poisson(p*trialCount.Rate);
		}

		public static Poisson TrialCountAverageConditional(int sample, double p, Poisson trialCount)
		{
			if (trialCount.Precision != 1) throw new NotImplementedException("precision != 1 not implemented");
			// p(n) =propto lambda^n/n! exp(-lambda) nchoosek(n,x) p^x (1-p)^(n-x)
			//      =propto lambda^n (1-p)^n 1/(n-x)!
			//      = Po(n-x; lambda*(1-p))
			// E[n] = x + lambda*(1-p)
			return new Poisson(sample + trialCount.Rate * (1-p)) / trialCount;
		}

		public static Poisson TrialCountAverageConditional(Poisson sample, double p)
		{
			if (sample.Precision != 0) throw new NotImplementedException("precision != 0 not implemented");
			// sum(x<=n) r^x nchoosek(n,x) p^x (1-p)^(n-x)
			// =propto (1-p)^n sum(x<=n) r^x nchoosek(n,x) (p/(1-p))^x
			// =propto (1-p)^n (1 + r p/(1-p))^n
			// =propto (1-p + r p)^n
			return new Poisson(1-p + sample.Rate*p, 0);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(trialCount) p(trialCount) factor(sample,trialCount,p)]/p(sample)</c>.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(double p, Discrete trialCount, Discrete result)
		{
			if (trialCount.IsPointMass) return SampleAverageConditional(p, trialCount.Point, result);
			// result must range from 0 to nMax
			if (result.Dimension < trialCount.Dimension) throw new ArgumentException("result.Dimension ("+result.Dimension+") < n.Dimension ("+trialCount.Dimension+")");
			Vector probs = result.GetWorkspace();
			double logp = Math.Log(p);
			double log1minusp = Math.Log(1-p);
			// p(sample=k) = sum_(n>=k) p(n) nchoosek(n,k) p^k (1-p)^(n-k)
			//             = (p/(1-p))^k 1/k! sum_(n>=k) p(n) n!/(n-k)! (1-p)^n
			for (int k = 0; k < result.Dimension; k++)
			{
				double s = 0.0;
				for (int n = k; n < trialCount.Dimension; n++)
				{
					s += trialCount[n] * Math.Exp(MMath.ChooseLn(n, k) + n*log1minusp);
				}
				probs[k] = Math.Exp(k*(logp-log1minusp))*s;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(p,trialCount) p(p,trialCount) factor(sample,trialCount,p)]/p(sample)</c>.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(Beta p, Discrete trialCount, Discrete result)
		{
			if (p.IsPointMass) return SampleAverageConditional(p.Point, trialCount, result);
			if (trialCount.IsPointMass) return SampleAverageConditional(p, trialCount.Point, result);
			// result must range from 0 to nMax
			if (result.Dimension < trialCount.Dimension) throw new ArgumentException("result.Dimension ("+result.Dimension+") < n.Dimension ("+trialCount.Dimension+")");
			Vector probs = result.GetWorkspace();
			double a = p.TrueCount;
			double b = p.FalseCount;
			// p(sample=k) = sum_(n>=k) p(n) nchoosek(n,k) p^k (1-p)^(n-k)
			//             = (p/(1-p))^k 1/k! sum_(n>=k) p(n) n!/(n-k)! (1-p)^n
			for (int k = 0; k < result.Dimension; k++)
			{
				double s = 0.0;
				for (int n = k; n < trialCount.Dimension; n++)
				{
					s += trialCount[n] * Math.Exp(MMath.ChooseLn(n, k) + MMath.GammaLn(b+n-k) - MMath.GammaLn(a+b+n));
				}
				probs[k] = Math.Exp(MMath.GammaLn(a+k))*s;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'trialCount'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'trialCount' as the random arguments are varied.
		/// The formula is <c>proj[p(trialCount) sum_(sample) p(sample) factor(sample,trialCount,p)]/p(trialCount)</c>.
		/// </para></remarks>
		public static Discrete TrialCountAverageConditional(Discrete sample, double p, Discrete result)
		{
			if (sample.IsPointMass) return TrialCountAverageConditional(sample.Point, p, result);
			// n must range from 0 to sampleMax
			if (result.Dimension < sample.Dimension) throw new ArgumentException("result.Dimension ("+result.Dimension+") < sample.Dimension ("+sample.Dimension+")");
			Vector probs = result.GetWorkspace();
			double logp = Math.Log(p);
			double log1minusp = Math.Log(1-p);
			// p(n) = sum_(k<=n) p(k) nchoosek(n,k) p^k (1-p)^(n-k)
			for (int n = 0; n < result.Dimension; n++)
			{
				double s = 0.0;
				for (int k = 0; k <= n; k++)
				{
					s += sample[k] * Math.Exp(MMath.ChooseLn(n, k) + k*(logp - log1minusp));
				}
				probs[n] = Math.Exp(n*log1minusp)*s;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'trialCount'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'trialCount' as the random arguments are varied.
		/// The formula is <c>proj[p(trialCount) sum_(sample,p) p(sample,p) factor(sample,trialCount,p)]/p(trialCount)</c>.
		/// </para></remarks>
		public static Discrete TrialCountAverageConditional(Discrete sample, Beta p, Discrete result)
		{
			if (p.IsPointMass) return TrialCountAverageConditional(sample, p.Point, result);
			if (sample.IsPointMass) return TrialCountAverageConditional(sample.Point, p, result);
			// n must range from 0 to sampleMax
			if (result.Dimension < sample.Dimension) throw new ArgumentException("result.Dimension ("+result.Dimension+") < sample.Dimension ("+sample.Dimension+")");
			Vector probs = result.GetWorkspace();
			double a = p.TrueCount;
			double b = p.FalseCount;
			// p(n) = sum_(k<=n) p(k) nchoosek(n,k) p^k (1-p)^(n-k)
			for (int n = 0; n < result.Dimension; n++)
			{
				double s = 0.0;
				for (int k = 0; k <= n; k++)
				{
					s += sample[k] * Math.Exp(MMath.ChooseLn(n, k) + MMath.GammaLn(a+k) + MMath.GammaLn(b+n-k));
				}
				probs[n] = Math.Exp(-MMath.GammaLn(a+b+n))*s;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'trialCount'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'trialCount' conditioned on the given values.
		/// </para></remarks>
		public static Discrete TrialCountAverageConditional(int sample, double p, Discrete result)
		{
			Vector probs = result.GetWorkspace();
			double logp = Math.Log(p);
			double log1minusp = Math.Log(1-p);
			// p(n) = nchoosek(n,k) p^k (1-p)^(n-k)
			for (int n = 0; n < result.Dimension; n++)
			{
				probs[n] = Math.Exp(MMath.ChooseLn(n, sample) + sample*logp + (n-sample)*log1minusp);
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'trialCount'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'trialCount' as the random arguments are varied.
		/// The formula is <c>proj[p(trialCount) sum_(p) p(p) factor(sample,trialCount,p)]/p(trialCount)</c>.
		/// </para></remarks>
		public static Discrete TrialCountAverageConditional(int sample, Beta p, Discrete result)
		{
			if (p.IsPointMass) return TrialCountAverageConditional(sample, p.Point, result);
			Vector probs = result.GetWorkspace();
			double a = p.TrueCount;
			double b = p.FalseCount;
			// p(n) = nchoosek(n,k) p^k (1-p)^(n-k)
			for (int n = 0; n < result.Dimension; n++)
			{
				probs[n] = Math.Exp(MMath.ChooseLn(n, sample) + MMath.GammaLn(b+n-sample) - MMath.GammaLn(a+b+n));
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'p'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>The outgoing EP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'p' conditioned on the given values.
		/// </para></remarks>
		public static Beta PAverageConditional(int sample, int trialCount)
		{
			// sample ranges from 0 to trialCount.
			return new Beta(sample+1, trialCount-sample+1);
		}
		public static Beta PAverageConditional(Discrete sample, int trialCount)
		{
			if (sample.IsPointMass) return PAverageConditional(sample.Point, trialCount);
			throw new NotSupportedException("Random sample is not yet implemented");
		}


		//-- VMP -----------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Constant value for 'p'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,trialCount,p))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(int sample, double p, int trialCount)
		{
			return LogAverageFactor(sample, p, trialCount);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(p) p(p) log(factor(sample,trialCount,p))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(int sample, Beta p, int trialCount)
		{
			double eLogP, eLogOneMinusP;
			p.GetMeanLogs(out eLogP, out eLogOneMinusP);
			return MMath.ChooseLn(trialCount, sample) + sample*eLogP + (trialCount-sample)*eLogOneMinusP;
		}

		/// <summary>
		/// VMP message to 'p'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>The outgoing VMP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'p' conditioned on the given values.
		/// </para></remarks>
		public static Beta PAverageLogarithm(int sample, int trialCount)
		{
			return PAverageConditional(sample, trialCount);
		}

		/// <summary>
		/// VMP message to 'p'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <returns>The outgoing VMP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'p'.
		/// The formula is <c>exp(sum_(trialCount) p(trialCount) log(factor(sample,trialCount,p)))</c>.
		/// </para></remarks>
		public static Beta PAverageLogarithm(int sample, Discrete trialCount)
		{
			// log f(sample,n,p) = sample*log(p) + (n-sample)*log(1-p) + const.
			return new Beta(sample+1, trialCount.GetMean()-sample+1);
		}

		/// <summary>
		/// VMP message to 'p'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>The outgoing VMP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'p'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,trialCount,p)))</c>.
		/// </para></remarks>
		public static Beta PAverageLogarithm(Discrete sample, int trialCount)
		{
			// log f(sample,n,p) = sample*log(p) + (n-sample)*log(1-p) + const.
			double mSample = sample.GetMean();
			return new Beta(mSample+1, trialCount-mSample+1);
		}

		/// <summary>
		/// VMP message to 'p'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Incoming message from 'trialCount'.</param>
		/// <returns>The outgoing VMP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'p'.
		/// The formula is <c>exp(sum_(sample,trialCount) p(sample,trialCount) log(factor(sample,trialCount,p)))</c>.
		/// </para></remarks>
		public static Beta PAverageLogarithm(Discrete sample, Discrete trialCount)
		{
			// log f(sample,n,p) = sample*log(p) + (n-sample)*log(1-p) + const.
			double mSample = sample.GetMean();
			return new Beta(mSample+1, trialCount.GetMean()-mSample+1);
		}
	}
}
