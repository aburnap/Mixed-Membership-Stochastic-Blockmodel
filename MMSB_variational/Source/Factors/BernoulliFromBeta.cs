// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Bernoulli"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Bernoulli), "Sample", typeof(double))]
	[FactorMethod(new String[] { "Sample", "ProbTrue" }, typeof(Factor), "Bernoulli")]
	[Quality(QualityBand.Mature)]
	public static class BernoulliFromBetaOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,probTrue))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool sample, double probTrue)
		{
			return sample ? Math.Log(probTrue) : Math.Log(1 - probTrue);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,probTrue))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli sample, [Fresh] Bernoulli to_sample)
		{
			return sample.GetLogAverageOf(to_sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probTrue">Incoming message from 'probTrue'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(probTrue) p(probTrue) factor(sample,probTrue))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool sample, Beta probTrue)
		{
			Bernoulli to_sample = SampleAverageConditional(probTrue);
			return to_sample.GetLogProb(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,probTrue) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli sample, double probTrue)
		{
			return 0.0;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probTrue">Incoming message from 'probTrue'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(probTrue) p(probTrue) factor(sample,probTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool sample, Beta probTrue)
		{
			return LogAverageFactor(sample, probTrue);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,probTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool sample, double probTrue)
		{
			return LogAverageFactor(sample, probTrue);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="probTrue">Incoming message from 'probTrue'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,probTrue) p(sample,probTrue) factor(sample,probTrue) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli sample, Beta probTrue)
		{
			return 0.0;
		}

		/// <summary>
		/// Gibbs message to 'sample'
		/// </summary>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing Gibbs message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli SampleConditional(double probTrue)
		{
			return new Bernoulli(probTrue);
		}
		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli SampleAverageConditional(double probTrue)
		{
			return SampleConditional(probTrue);
		}
		/// <summary>
		/// Gibbs message to 'probTrue'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <returns>The outgoing Gibbs message to the 'probTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'probTrue' conditioned on the given values.
		/// </para></remarks>
		public static Beta ProbTrueConditional(bool sample)
		{
			if (sample) {
				return new Beta(2, 1);
			} else {
				return new Beta(1, 2);
			}
		}
		/// <summary>
		/// EP message to 'probTrue'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <returns>The outgoing EP message to the 'probTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'probTrue' conditioned on the given values.
		/// </para></remarks>
		public static Beta ProbTrueAverageConditional(bool sample)
		{
			return ProbTrueConditional(sample);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="probTrue">Incoming message from 'probTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(probTrue) p(probTrue) factor(sample,probTrue)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probTrue"/> is not a proper distribution</exception>
		public static Bernoulli SampleAverageConditional([SkipIfUniform] Beta probTrue)
		{
			if (probTrue.IsPointMass) {
				return new Bernoulli(probTrue.Point);
			} else if (!probTrue.IsProper()) throw new ImproperMessageException(probTrue);
			else {
				// p(x=true) = trueCount/total
				// p(x=false) = falseCount/total
				// log(p(x=true)/p(x=false)) = log(trueCount/falseCount)
				return Bernoulli.FromLogOdds(Math.Log(probTrue.TrueCount / probTrue.FalseCount));
			}
		}

		/// <summary>
		/// EP message to 'probTrue'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="probTrue">Incoming message from 'probTrue'.</param>
		/// <returns>The outgoing EP message to the 'probTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'probTrue' as the random arguments are varied.
		/// The formula is <c>proj[p(probTrue) sum_(sample) p(sample) factor(sample,probTrue)]/p(probTrue)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Beta ProbTrueAverageConditional([SkipIfUniform] Bernoulli sample, Beta probTrue)
		{
			// this code is similar to DiscreteFromDirichletOp.PAverageConditional()
			if (probTrue.IsPointMass) {
				return Beta.Uniform();
			}
			if (sample.IsPointMass) {
				// shortcut
				return ProbTrueConditional(sample.Point);
			}
			if (!probTrue.IsProper()) throw new ImproperMessageException(probTrue);
			// q(x) is the distribution stored in this.X.
			// q(p) is the distribution stored in this.P.
			// f(x,p) is the factor.
			// Z = sum_x q(x) int_p f(x,p)*q(p) = q(false)*E[1-p] + q(true)*E[p]
			// Ef[p] = 1/Z sum_x q(x) int_p p*f(x,p)*q(p) = 1/Z (q(false)*E[p(1-p)] + q(true)*E[p^2])
			// Ef[p^2] = 1/Z sum_x q(x) int_p p^2*f(x,p)*q(p) = 1/Z (q(false)*E[p^2(1-p)] + q(true)*E[p^3])
			// var_f(p) = Ef[p^2] - Ef[p]^2
			double mo = probTrue.GetMean();
			double m2o = probTrue.GetMeanSquare();
			double pT = sample.GetProbTrue();
			double pF = sample.GetProbFalse();
			double Z = pF * (1 - mo) + pT * mo;
			double m = pF * (mo - m2o) + pT * m2o;
			m = m / Z;
			if (!Beta.AllowImproperSum) {
				if (pT < 0.5) {
					double inc = probTrue.TotalCount * (mo / m - 1);
					return new Beta(1, 1 + inc);
				} else {
					double inc = probTrue.TotalCount * ((1 - mo) / (1 - m) - 1);
					return new Beta(1 + inc, 1);
				}
			} else {
				double m3o = probTrue.GetMeanCube();
				double m2 = pF * (m2o - m3o) + pT * m3o;
				m2 = m2 / Z;
				Beta result = Beta.FromMeanAndVariance(m, m2 - m * m);
				result.SetToRatio(result, probTrue);
				return result;
			}
		}

		//-- VMP -------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="probTrue">Incoming message from 'probTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,probTrue) p(sample,probTrue) log(factor(sample,probTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="probTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] Bernoulli sample, [Proper] Beta probTrue)
		{
			if (sample.IsPointMass) return AverageLogFactor(sample.Point, probTrue);
			double eLogP, eLog1MinusP;
			probTrue.GetMeanLogs(out eLogP, out eLog1MinusP);
			double p = sample.GetProbTrue();
			return p * eLogP + (1 - p) * eLog1MinusP;
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,probTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] Bernoulli sample, double probTrue)
		{
			if (sample.IsPointMass) return AverageLogFactor(sample.Point, probTrue);
			return AverageLogFactor(sample, Beta.PointMass(probTrue));
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probTrue">Incoming message from 'probTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(probTrue) p(probTrue) log(factor(sample,probTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor(bool sample, [Proper] Beta probTrue)
		{
			double eLogP, eLog1MinusP;
			probTrue.GetMeanLogs(out eLogP, out eLog1MinusP);
			return sample ? eLogP : eLog1MinusP;
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,probTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool sample, double probTrue)
		{
			return sample ? Math.Log(probTrue) : Math.Log(1 - probTrue);
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli SampleAverageLogarithm(double probTrue)
		{
			return SampleConditional(probTrue);
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="probTrue">Incoming message from 'probTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(probTrue) p(probTrue) log(factor(sample,probTrue)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probTrue"/> is not a proper distribution</exception>
		public static Bernoulli SampleAverageLogarithm([SkipIfUniform] Beta probTrue)
		{
			if (probTrue.IsPointMass) {
				return new Bernoulli(probTrue.Point);
			} else if (!probTrue.IsProper()) throw new ImproperMessageException(probTrue);
			else {
				// E[x*log(p) + (1-x)*log(1-p)] = x*E[log(p)] + (1-x)*E[log(1-p)]
				// p(x=true) = exp(E[log(p)])/(exp(E[log(p)]) + exp(E[log(1-p)]))
				// log(p(x=true)/p(x=false)) = E[log(p)] - E[log(1-p)] = digamma(trueCount) - digamma(falseCount)
				return Bernoulli.FromLogOdds(MMath.Digamma(probTrue.TrueCount) - MMath.Digamma(probTrue.FalseCount));
			}
		}
		/// <summary>
		/// VMP message to 'probTrue'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <returns>The outgoing VMP message to the 'probTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'probTrue' conditioned on the given values.
		/// </para></remarks>
		public static Beta ProbTrueAverageLogarithm(bool sample)
		{
			return ProbTrueConditional(sample);
		}
		/// <summary>
		/// VMP message to 'probTrue'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'probTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'probTrue'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,probTrue)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Beta ProbTrueAverageLogarithm(Bernoulli sample)
		{
			// E[x*log(p) + (1-x)*log(1-p)] = E[x]*log(p) + (1-E[x])*log(1-p)
			double ex = sample.GetProbTrue();
			return new Beta(1 + ex, 2 - ex);
		}
	}
}
