// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.BernoulliFromLogOdds"/>, given random arguments to the function.
	/// </summary> 
	[FactorMethod(typeof(Factor), "BernoulliFromLogOdds", Default = true)]
	[Quality(QualityBand.Experimental)]
	[Obsolete("Use LogisticOp_JJ96")]
	public static class BernoulliFromLogOddsOp_JJ96
	{
		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Fixed value for sample</param>
		/// <param name="logOdds">Incoming message from logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double AverageLogFactor(bool sample, [Proper, SkipIfUniform] Gaussian logOdds)
		{
			if (logOdds.IsUniform()) return 0.0;
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			double t = Math.Sqrt(m * m + v);
			double s = sample ? 1 : -1;
			return MMath.LogisticLn(t) + (s * m - t) / 2;
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Incoming message from sample</param>
		/// <param name="logOdds">Incoming message from logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double AverageLogFactor(Bernoulli sample, [Proper, SkipIfUniform] Gaussian logOdds)
		{
			if (logOdds.IsUniform()) return 0.0;
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			double t = Math.Sqrt(m * m + v);
			double s = 2 * sample.GetProbTrue() - 1;  // probTrue - probFalse
			return MMath.LogisticLn(t) + (s * m - t) / 2;
		}

		/// <summary>
		/// VMP message to LogOdds
		/// </summary>
		/// <param name="sample">Fixed value for sample</param>
		/// <param name="logOdds">Incoming message from logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'logOdds'.
		/// The formula is <c>int log(f(logOdds,x)) q(x) dx</c> where <c>x = (sample)</c>.
		/// </para></remarks>
		public static Gaussian LogOddsAverageLogarithm(bool sample, [Proper, SkipIfUniform] Gaussian logOdds)
		{
			// sigma(x) >= sigma(t) exp((x-t)/2 - a/2*(x^2 - t^2))
			if (logOdds.IsUniform()) return logOdds;
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			double t = Math.Sqrt(m * m + v);
			double lambda = (t == 0) ? 0.25 : Math.Tanh(t / 2) / (2 * t);
			return Gaussian.FromNatural(sample ? 0.5 : -0.5, lambda);
		}

		/// <summary>
		/// VMP message to LogOdds
		/// </summary>
		/// <param name="sample">Incoming message from sample</param>
		/// <param name="logOdds">Incoming message from logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'logOdds'.
		/// The formula is <c>int log(f(logOdds,x)) q(x) dx</c> where <c>x = (sample)</c>.
		/// </para></remarks>
		public static Gaussian LogOddsAverageLogarithm(Bernoulli sample, [Proper, SkipIfUniform] Gaussian logOdds)
		{
			if (logOdds.IsUniform()) return logOdds;
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			double t = Math.Sqrt(m * m + v);
			double lambda = (t == 0) ? 0.25 : Math.Tanh(t / 2) / (2 * t);
			return Gaussian.FromNatural(sample.GetProbTrue() - 0.5, lambda);
		}
	}

    /// <summary>
    /// Provides outgoing messages for <see cref="Factor.BernoulliFromLogOdds"/>, given random arguments to the function.
    /// Uses a generalisation of the bound from Saul and Jordan (1999). 
    /// </summary>
	[FactorMethod(typeof(Factor), "BernoulliFromLogOdds")]
	[Quality(QualityBand.Experimental)]
	[Obsolete("Use LogisticOp_SJ99")]
	public static class BernoulliFromLogOddsOp_SJ99
	{
        /// <summary>
		/// VMP message to LogOdds
		/// </summary>
		/// <param name="sample">Incoming message from sample</param>
		/// <param name="logOdds">Incoming message from logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'logOdds'.
		/// The formula is <c>int log(f(logOdds,x)) q(x) dx</c> where <c>x = (sample)</c>.
		/// </para></remarks>
		public static Gaussian LogOddsAverageLogarithm(bool sample, Gaussian logOdds)
		{
			// This is the non-conjugate VMP update using the Saul and Jordan (1999) bound.
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			double a = 0.5;
			// TODO: use a buffer to store the value of 'a', so it doesn't need to be re-optimised each time.
			for (int iter = 0; iter < 10; iter++) {
				double aOld = a;
				a = MMath.Logistic(m + (1-2*a)*v*0.5);
				if (Math.Abs(a - aOld) < 1e-8) break;
			}
			double sa = MMath.Logistic(m + (1-2*a)*v*0.5);
			double vf = 1/(a*a + (1-2*a)*sa);
			double mf = m + vf*(sample ? 1-sa : sa);
			return Gaussian.FromMeanAndVariance(mf, vf);
		}

        public static double AverageLogFactor(Bernoulli sample, Gaussian logOdds)
        {
            // This is the non-conjugate VMP update using the Saul and Jordan (1999) bound.
            double m, v;
            logOdds.GetMeanAndVariance(out m, out v);
            double a = 0.5;
            // TODO: use a buffer to store the value of 'a', so it doesn't need to be re-optimised each time.
            for (int iter = 0; iter < 10; iter++)
            {
                double aOld = a;
                a = MMath.Logistic(m + (1 - 2 * a) * v * 0.5);
                if (Math.Abs(a - aOld) < 1e-8) break;
            }
            return sample.GetProbTrue() * m - .5 * a * a * v - MMath.Log1PlusExp(m + (1 - 2 * a) * v * 0.5);
        }


	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.BernoulliFromLogOdds"/>, given random arguments to the function.
	/// Performs KL minimisation using gradient matching, a distributed gradient descent algorithm. 
	/// </summary>
	[FactorMethod(typeof(Factor), "BernoulliFromLogOdds")]
	[Quality(QualityBand.Experimental)]
	[Obsolete("Use LogisticOp")]
	public static class BernoulliFromLogOddsOp
	{
		/// <summary>
		/// EP message to 'logOdds'.
		/// </summary>
		/// <param name="sample">Constant value for sample.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'logOdds' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the moment matched Gaussian approximation to the factor.
		/// </para></remarks>
		public static Gaussian LogOddsAverageConditional(bool sample, [SkipIfUniform] Gaussian logOdds)
		{
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			double s = sample ? 1 : -1;
			m *= s;
			if (m + 1.5*v < -38) {
				double beta2 = Math.Exp(m + 1.5*v);
				return Gaussian.FromMeanAndVariance(s*(m + v), v*(1 - v*beta2)) / logOdds;
			}
			double sigma0 = MMath.LogisticGaussian(m, v);
			double sigma1 = MMath.LogisticGaussianDerivative(m, v);
			double sigma2 = MMath.LogisticGaussianDerivative2(m, v);
			double alpha, beta;
			alpha = sigma1 / sigma0;
			if (Double.IsNaN(alpha)) throw new Exception("alpha is NaN");
			if (m + 2*v < -19) {
				beta = Math.Exp(3*m + 2.5*v)/(sigma0*sigma0);
			} else {
				//beta = (sigma1*sigma1 - sigma2*sigma0)/(sigma0*sigma0);
				beta = alpha*alpha - sigma2/sigma0;
			}
			if (Double.IsNaN(beta)) throw new Exception("beta is NaN");
			double m2 = s * (m + v * alpha);
			double v2 = v*(1 - v*beta);
			if (v2 > v) throw new Exception("v2 > v");
			if (v2 < 0) throw new Exception("v2 < 0");
			return Gaussian.FromMeanAndVariance(m2, v2) / logOdds;
		}

		/// <summary>
		/// EP message to 'logOdds'.
		/// </summary>
		/// <param name="sample">Incoming message from sample.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'logOdds' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the moment matched Gaussian approximation to the factor.
		/// </para></remarks>
		public static Gaussian LogOddsAverageConditional(Bernoulli sample, [SkipIfUniform] Gaussian logOdds)
		{
			Gaussian toLogOddsT = LogOddsAverageConditional(true, logOdds);
			double logWeightT = LogAverageFactor(true, logOdds) + sample.GetLogProbTrue();
			Gaussian toLogOddsF = LogOddsAverageConditional(false, logOdds);
			double logWeightF = LogAverageFactor(false, logOdds) + sample.GetLogProbFalse();
			double maxWeight = Math.Max(logWeightT, logWeightF);
			logWeightT -= maxWeight;
			logWeightF -= maxWeight;
			Gaussian result = new Gaussian();
			result.SetToSum(Math.Exp(logWeightT), toLogOddsT * logOdds, Math.Exp(logWeightF), toLogOddsF * logOdds);
			result /= logOdds;
			return result;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Constant value for 'logOdds'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool sample, double logOdds)
		{
			return MMath.LogisticLn(sample ? logOdds : -logOdds);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool sample, Gaussian logOdds)
		{
			return LogAverageFactor(sample, logOdds);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Constant value for 'logOdds'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool sample, double logOdds)
		{
			return LogAverageFactor(sample, logOdds);
		}

        /// <summary>
        /// Evidence message for EP.
        /// </summary>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli sample) { return 0.0; }

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool sample, [Proper] Gaussian logOdds)
		{
			if (logOdds.IsPointMass) return LogAverageFactor(sample, logOdds.Point);
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			return Math.Log(MMath.LogisticGaussian(sample ? m : -m, v));
		}

		/// <summary>
		/// EP message to 'sample'.
		/// </summary>
		/// <param name="logOdds">Incoming message from 'logOdds'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the moment matched Gaussian approximation to the factor.
		/// </para></remarks>
		public static Bernoulli SampleAverageConditional([Proper] Gaussian logOdds)
		{
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			return new Bernoulli(MMath.LogisticGaussian(m, v));
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Fixed value for sample</param>
		/// <param name="logOdds">Fixed value for logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double AverageLogFactor(bool sample, double logOdds)
		{
			if (sample) return MMath.LogisticLn(logOdds);
			else return MMath.LogisticLn(-logOdds);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Incoming message from sample</param>
		/// <param name="logOdds">Fixed value for logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (sample,logOdds)</c>.
		/// </para></remarks>
		public static double AverageLogFactor(Bernoulli sample, double logOdds)
		{
			if (sample.IsPointMass) return AverageLogFactor(sample.Point, logOdds);
			// probTrue*log(sigma(logOdds)) + probFalse*log(sigma(-logOdds))
			// = -log(1+exp(-logOdds)) + probFalse*(-logOdds)
			// = probTrue*logOdds - log(1+exp(logOdds))
			if (logOdds >= 0) {
				double probFalse = sample.GetProbFalse();
				return -probFalse * logOdds - MMath.Log1PlusExp(-logOdds);
			} else {
				double probTrue = sample.GetProbTrue();
				return probTrue * logOdds - MMath.Log1PlusExp(logOdds);
			}
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// </para></remarks>
		public static double AverageLogFactor(bool sample, [Proper, SkipIfUniform] Gaussian logOdds)
		{
			// f(sample,logOdds) = exp(sample*logOdds)/(1 + exp(logOdds))
			// log f(sample,logOdds) = sample*logOdds - log(1 + exp(logOdds))
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			return (sample ? 1.0 : 0.0) * m - MMath.Log1PlusExpGaussian(m, v);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// </para></remarks>
		public static double AverageLogFactor(Bernoulli sample, [Proper, SkipIfUniform] Gaussian logOdds)
		{
			// f(sample,logOdds) = exp(sample*logOdds)/(1 + exp(logOdds))
			// log f(sample,logOdds) = sample*logOdds - log(1 + exp(logOdds))
			double m, v;
			logOdds.GetMeanAndVariance(out m, out v);
			return sample.GetProbTrue() * m - MMath.Log1PlusExpGaussian(m, v);
		}

		// Calculate \int f(x) dx over the whole real line. Uses a change of variable
		// x=cot(t) which maps the real line onto [0,pi] and use trapezium rule
		private static double ClenshawCurtisQuadrature(Converter<double, double> f, int CCFactor, int numIntervals)
		{
			double intervalWidth = Math.PI / (double)numIntervals;
			double sum = 0;
			for (double x = intervalWidth; x < Math.PI; x += intervalWidth) {
				double sinX = Math.Sin(x);
				sum += f(CCFactor / Math.Tan(x)) / (sinX * sinX);
			}
			return CCFactor * intervalWidth * sum;
		}

		/// <summary>
		/// VMP message to sample
		/// </summary>
		/// <param name="logOdds">Incoming message from logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (logOdds)</c>.
		/// </para></remarks>
		public static Bernoulli SampleAverageLogarithm([Proper] Gaussian logOdds)
		{
			return Bernoulli.FromLogOdds(logOdds.GetMean());
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="logOdds">Fixed value for logOdds</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (logOdds)</c>.
		/// </para></remarks>
		public static Bernoulli SampleAverageLogarithm(double logOdds)
		{
			return Bernoulli.FromLogOdds(logOdds);
		}

		/// <summary>
		/// Gradient matching VMP message from factor to logOdds variable
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Incoming message. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_LogOdds">Previous message sent, used for damping</param>
		/// <returns>The outgoing VMP message.</returns>
		/// <remarks><para>
		/// The outgoing message is the Gaussian approximation to the factor which results in the 
		/// same derivatives of the KL(q||p) divergence with respect to the parameters of the posterior
		/// as if the true factor had been used.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="logOdds"/> is not a proper distribution</exception>
		public static Gaussian LogOddsAverageLogarithm(bool sample, [Proper, SkipIfUniform] Gaussian logOdds, Gaussian to_LogOdds)
		{
			double m, v; // prior mean and variance
			double s = sample ? 1 : -1;
			logOdds.GetMeanAndVariance(out m, out v);
			// E = \int q log f dx
			// Match gradients
			double dEbydm = s * MMath.LogisticGaussian(-s * m, v);
			double dEbydv = -.5 * MMath.LogisticGaussianDerivative(s * m, v);
			double prec = -2.0 * dEbydv;
			double meanTimesPrec = m * prec + dEbydm;
			Gaussian result = Gaussian.FromNatural(meanTimesPrec, prec);
			double step = Rand.Double() * 0.5; // random damping helps convergence, especially with parallel updates
			if (step != 1.0) {
				result.Precision = step * result.Precision + (1 - step) * to_LogOdds.Precision;
				result.MeanTimesPrecision = step * result.MeanTimesPrecision + (1 - step) * to_LogOdds.MeanTimesPrecision;
			}
			return result;
		}

		/// <summary>
		/// Gradient matching VMP message from factor to logOdds variable
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="logOdds">Incoming message. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Previous message sent, used for damping</param>
		/// <returns>The outgoing VMP message.</returns>
		/// <remarks><para>
		/// The outgoing message is the Gaussian approximation to the factor which results in the 
		/// same derivatives of the KL(q||p) divergence with respect to the parameters of the posterior
		/// as if the true factor had been used.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="logOdds"/> is not a proper distribution</exception>
		public static Gaussian LogOddsAverageLogarithm(Bernoulli sample, [Proper, SkipIfUniform] Gaussian logOdds, Gaussian result)
		{
			if (sample.IsUniform()) return Gaussian.Uniform();
			throw new NotImplementedException("BernoulliFromLogOdds with non-observed output is not yet implemented");
		}
	}
}
