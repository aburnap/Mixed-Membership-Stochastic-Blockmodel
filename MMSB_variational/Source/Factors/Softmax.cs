// (C) Copyright 2009 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Collections;

using GaussianArray = MicrosoftResearch.Infer.Distributions.DistributionStructArray<MicrosoftResearch.Infer.Distributions.Gaussian, double>;



namespace MicrosoftResearch.Infer.Factors
{

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the simple first order Taylor series expansion from Blei et al. 06, followed by
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Experimental)]
	[Buffers("A")]
	public static class BouchardSoftmaxOp_Sparse
	{
		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static double AInit()
		{
			return .5;
		}

		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static double A<GaussianList>([SkipIfAnyUniform] GaussianList x, double a) where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			double old_a = a + 1;
			int iter = 0;
			//var acheck = BouchardSoftmaxOp.A(x, a);
			while (Math.Abs(a - old_a) > 1e-2 && iter < 100) {
				old_a = a;
				Func<Gaussian, double[]> f = xi => {
					double m, v;
					xi.GetMeanAndVariance(out m, out v);
					double aMinusM = a - m;
					double t = Math.Sqrt(aMinusM * aMinusM + v);
					double grad_t = aMinusM / t;
					double logisticT = MMath.Logistic(t);
					double hess_t = 1.0 / t - aMinusM * aMinusM / (t * t * t);
					return new double[] { -(1.0 + grad_t) / 2.0 + grad_t * logisticT, // grad
                            hess_t * (logisticT - .5) + grad_t * grad_t * logisticT * MMath.Logistic(-t) };
				};
				// todo when we go to .NET 4.0: use tuple instead of 2-array
				double[] gradHess = f.Map(x).EnumerableReduce(new double[] { 0, 0 },
								(res, current) => { res[0] += current[0]; res[1] += current[1]; return res; },
								(res, commonValue, commonValueCount) => {
									res[0] += commonValue[0] * commonValueCount;
									res[1] += commonValue[1] * commonValueCount;
									return res;
								});
				double grad = gradHess[0] + 1.0;
				double hess = gradHess[1];
				if (hess < .1)
					hess = .1;
				a -= grad / hess;
				iter++;
			}

			return a;
		}

		public static double logSumExp<GaussianList>(GaussianList x, double a)
				where GaussianList : IList<Gaussian>
		{
			return a + x.EnumerableSum(i => {
				double m, v;
				i.GetMeanAndVariance(out m, out v);
				double t = Math.Sqrt((m - a) * (m - a) + v);
				return (m - a - t) / 2.0 - MMath.LogisticLn(-t);
			});
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax, double a)
				where GaussianList : IList<Gaussian>
		{
			int K = softmax.Dimension;
			Func<Gaussian, double, double> f = (i, j) => i.GetMean() * (j - 1.0);
			double innerProduct = f.Map(x, softmax.PseudoCount).EnumerableSum(i => i);
			double sum_n = softmax.PseudoCount.EnumerableSum(i => i - 1.0);
			return innerProduct - sum_n * logSumExp(x, a) - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax); ;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, double a, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			double sum = a + x.EnumerableSum(i => {
				double m, v;
				i.GetMeanAndVariance(out m, out v);
				double t = Math.Sqrt((m - a) * (m - a) + v);
				return (m - a - t) / 2.0 - MMath.LogisticLn(-t);
			});
			var meanLog = SparseVector.Zero(x.Count);
			Func<Gaussian, double> f = i => i.GetMean() - sum;
			meanLog.SetTo(f.Map(x));
			result.SetMeanLog(meanLog);
			return result;
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax,
				[SkipIfAllUniform] GaussianList x,
				double a,
				GaussianList result)
				where GaussianList : IList<Gaussian>
		{
			double step = 1.0;
			int K = x.Count;
			double total = softmax.PseudoCount.Sum() - K;
			Func<Gaussian, double, Gaussian> f = (i, j) => {
				double m, v;
				i.GetMeanAndVariance(out m, out v);
				double t = Math.Sqrt((m - a) * (m - a) + v);
				double lambda = (1.0 / (1.0 + Math.Exp(-t)) - .5) / (2.0 * t);
				return Gaussian.FromNatural(j - 1.0 - total * (.5 - 2.0 * a * lambda), 2.0 * total * lambda);
				//return Gaussian.FromNatural(j - 1.0 - (.5 - 2.0 * a * lambda), 2.0 * lambda);
			};
			if (step == 1.0) {
				result.SetTo(f.Map(x, softmax.PseudoCount));
			} else {
				Func<Gaussian, Gaussian, Gaussian> damp = (i, j) => (i ^ step) * (j ^ (1.0 - step));
				result.SetTo(damp.Map(f.Map(x, softmax.PseudoCount), result));
			}
			return result;
		}

	}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the simple first order Taylor series expansion from Blei et al. 06, followed by
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Experimental)]
	[Buffers("A")]
	public static class BouchardSoftmaxOp
	{
		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static double AInit()
		{
			return .5;
		}

		static double[] GetT(double a, double[] m, double[] v)
		{
			var result = new double[m.Length];
			for (int i = 0; i < m.Length; i++) {
				result[i] = Math.Sqrt((m[i] - a) * (m[i] - a) + v[i]);
			}
			return result;
		}

		static double[] GetLambda(double[] t)
		{
			var result = new double[t.Length];
			for (int i = 0; i < t.Length; i++) {
				result[i] = (1.0 / (1.0 + Math.Exp(-t[i])) - .5) / (2.0 * t[i]);
			}
			return result;
		}

		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static double A<GaussianList>([SkipIfAnyUniform] GaussianList x, double a) where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			double old_a = a + 1;
			int iter = 0;
			while (Math.Abs(a - old_a) > 1e-2 && iter < 100) {
				old_a = a;
				double grad = 1.0, hess = 0.0;
				for (int k = 0; k < K; k++) {
					double aMinusM = a - m[k];
					double t = Math.Sqrt(aMinusM * aMinusM + v[k]);
					double grad_t = aMinusM / t;
					double hess_t = 1.0 / t - aMinusM * aMinusM / (t * t * t);
					double logisticT = MMath.Logistic(t);
					grad += -(1.0 + grad_t) / 2.0 + grad_t * logisticT;
					hess += hess_t * (logisticT - .5) + grad_t * grad_t * logisticT * MMath.Logistic(-t);
				}
				if (hess < .1)
					hess = .1;
				a -= grad / hess;
				iter++;
			}
			return a;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax, double a)
				 where GaussianList : IList<Gaussian>
		{
			int K = softmax.Dimension;
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var t = GetT(a, m, v);
			double logSumExp = a;

			for (int k = 0; k < K; k++) {
				logSumExp += (m[k] - a - t[k]) / 2.0 - MMath.LogisticLn(-t[k]);
			}
			var ns = softmax.PseudoCount - 1.0;
			double innerProduct = ns.Inner(Vector.FromArray(m));
			double sum_n = ns.Sum();
			return innerProduct - sum_n * logSumExp - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax); ;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, double a, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var t = GetT(a, m, v);
			double sum = a;
			for (int k = 0; k < K; k++) {
				sum += (m[k] - a - t[k]) / 2.0 - MMath.LogisticLn(-t[k]);
			}
			var meanLog = Vector.Constant(K, -sum);
			meanLog += Vector.FromArray(m);

			var sanityCheck = MMath.Softmax(m);
			var sanityCheck2 = meanLog.Select(Math.Exp).ToArray();
			for (int k = 0; k < K; k++) {
				if (sanityCheck2[k] > sanityCheck[k])
					throw new ApplicationException("doh");
			}
			result.SetMeanLog(meanLog);
			return result;
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax,
				[SkipIfAllUniform] GaussianList x,
				double a,
				GaussianList result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			double total = softmax.TotalCount - K;
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var t = GetT(a, m, v);
			var lambda = GetLambda(t);

			for (int k = 0; k < K; k++) {
				var rk = result[k];
				rk.Precision = 2 * total * lambda[k];
				rk.MeanTimesPrecision = softmax.PseudoCount[k] - 1.0 - total * (.5 - 2.0 * a * lambda[k]);
				result[k] = rk;
			}
			return result;
		}

	}


	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the simple first order Taylor series expansion from Blei and Lafferty (2006), followed by
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	public static class SoftmaxOp_BL06_LBFGS
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			double[] ms, vs;
			GetMeanAndVariance(x, out ms, out vs);
			double logSumExp = LogSumExpMPlusHalfV(ms, vs);
			var ns = softmax.PseudoCount - 1.0;
			double sum_n = ns.Sum();
			return ns.Inner(Vector.FromArray(ms)) - sum_n * logSumExp - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax); ;
		}

		/// <summary>
		/// Helper function to get the means and variances of a list of Gaussians
		/// </summary>
		internal static void GetMeanAndVariance(IList<Gaussian> x, out double[] m, out double[] v)
		{
			int K = x.Count;
			m = new double[K];
			v = new double[K];
			for (int k = 0; k < K; k++) {
				double tmp_m, tmp_v;
				x[k].GetMeanAndVariance(out tmp_m, out tmp_v);
				m[k] = tmp_m;
				v[k] = tmp_v;
			}
		}

		/// <summary>
		/// Helper function to get the means and variances of a list of Gaussians
		/// </summary>
		internal static void GetMeanAndVariance(IList<Gaussian> x, out Vector m, out Vector v)
		{
			double tmp_m, tmp_v;
			int K = x.Count;
			bool isSparse = x.IsSparse();
			m = Vector.Zero(K);
			v = Vector.Zero(K);

			for (int k = 0; k < K; k++) {
				x[k].GetMeanAndVariance(out tmp_m, out tmp_v);
				m[k] = tmp_m;
				v[k] = tmp_v;
			}
		}

		/// <summary>
		/// Helper function to calculation log sum_k exp(m_k+v_k/2), which is an upper bound
		/// on E[log sum_k exp(x_k)] when x_k ~ N(m_k,v_k)
		/// </summary>
		private static double LogSumExpMPlusHalfV(double[] ms, double[] vs)
		{
			var mPlusHalfV = Vector.FromArray(ms);
			var halfV = Vector.FromArray(vs);
			halfV.Scale(.5);
			mPlusHalfV += halfV;
			return MMath.LogSumExp(mPlusHalfV);
		}

		/// <summary>
		/// Function to evaluate this factor's KL divergence contribution, and gradient (if grad is not null). 
		/// </summary>
		/// <param name="mu">Prior means</param>
		/// <param name="s2">Prior variances</param>
		/// <param name="x">x[1..K]: posterior mean, x[K+1..2K]: posterior log variance</param>
		/// <param name="ns">Dirichlet counts-1</param>
		/// <param name="grad">Vector to store the gradient in</param>
		/// <returns>KL divergence</returns>
		private static double GradientAndValueAtPoint(double[] mu, double[] s2, Vector x, Vector ns, Vector grad)
		{
			int K = x.Count / 2;
			var ms = new double[K];
			var vs = new double[K];
			for (int k = 0; k < K; k++) {
				ms[k] = x[k];
				vs[k] = Math.Exp(x[k + K]);
			}
			double sum_n = ns.Sum();
			double logSumExp = LogSumExpMPlusHalfV(ms, vs);
			double kl_value = sum_n * logSumExp;
			for (int k = 0; k < K; k++) {
				if (s2[k] > 0.0)
					kl_value += -.5 * (1.0 + x[K + k]) + .5 * (vs[k] + (ms[k] - mu[k]) * (ms[k] - mu[k])) / s2[k] + .5 * Math.Log(s2[k]) - ns[k] * ms[k];
				if (grad != null) {
					if (s2[k] == 0.0) {
						grad[k] = 0.0;
						grad[K + k] = 0.0;
					} else {
						double logistic_k = Math.Exp(ms[k] + vs[k] / 2 - logSumExp);
						grad[k] = (ms[k] - mu[k]) / s2[k] - ns[k] + sum_n * logistic_k;
						grad[K + k] = -.5 + .5 * vs[k] / s2[k] + .5 * vs[k] * sum_n * logistic_k;
					}
				}
			}

			return kl_value;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			// TM: this routine has been rewritten to allocate as little memory as possible.
			int K = x.Count;
			Vector meanLog = Vector.Zero(K);
			double[] meanPlusHalfV = new double[K];
			for (int i = 0; i < K; i++) {
				double m, v;
				x[i].GetMeanAndVariance(out m, out v);
				meanLog[i] = m;
				meanPlusHalfV[i] = m + 0.5 * v;
			}
			double sum = MMath.LogSumExp(meanPlusHalfV);
			for (int i = 0; i < K; i++) {
				meanLog[i] -= sum;
			}
			result.SetMeanLog(meanLog);
			return result;
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [Stochastic, SkipIfAllUniform] IList<Gaussian> x, GaussianList result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			var prior = new Gaussian[K];
			for (int k = 0; k < K; k++)
				prior[k] = x[k] / result[k];
			double[] mu, s2;
			GetMeanAndVariance(prior, out mu, out s2);
			double[] m, v;
			GetMeanAndVariance(x, out m, out v);
			var counts = softmax.PseudoCount - 1;
			int evalCounter = 0;
			// z[1..K]: posterior means, z[K+1..2K]: log posterior variances
			var z = Vector.Zero(K * 2);
			for (int k = 0; k < K; k++) {
				z[k] = m[k];
				z[K + k] = Math.Log(v[k]);
			}
			double startingValue = GradientAndValueAtPoint(mu, s2, z, counts, null);
			var s = new LBFGS(5);
			s.MaximumStep = 1e3;
			s.MaximumIterations = 100;
			s.Epsilon = 1e-10;
			s.convergenceCriteria = BFGS.ConvergenceCriteria.Objective;
			z = s.Run(z, 1.0, delegate(Vector y, ref Vector grad) { evalCounter++; return GradientAndValueAtPoint(mu, s2, y, counts, grad); });
			for (int k = 0; k < K; k++) {
				m[k] = z[k];
				v[k] = Math.Exp(z[K + k]);
				var rk = result[k];
				if (!prior[k].IsPointMass)
					rk.SetToRatio(Gaussian.FromMeanAndVariance(m[k], v[k]), prior[k]);
				result[k] = rk;
			}
			double endValue = GradientAndValueAtPoint(mu, s2, z, counts, null);
			//Console.WriteLine("Went from {0} to {1} in {2} steps, {3} evals", startingValue, endValue, s.IterationsPerformed, evalCounter);
			if (startingValue < endValue)
				Console.WriteLine("Warning: LBFGS resulted in an increased objective function");
			return result;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		//public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax)
		//{
		//    double logSumExp = LogSumExpMPlusHalfV(x);
		//    var ns = softmax.PseudoCount - 1.0;
		//    double innerProduct = x.ListZip(softmax.PseudoCount, (i, j) => i.GetMean() * (j - 1.0))
		//        .EnumerableSum(i => i);
		//    double sum_n = softmax.PseudoCount.EnumerableSum(i => i - 1.0);
		//    return innerProduct - sum_n * logSumExp - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax); ;
		//}

		///// <summary>
		///// Helper function to calculation log sum_k exp(m_k+v_k/2), which is an upper bound
		///// on E[log sum_k exp(x_k)] when x_k ~ N(m_k,v_k)
		///// </summary>
		//private static double LogSumExpMPlusHalfV(IList<Gaussian> x)
		//{
		//    return MMath.LogSumExpSparse(x.ListSelect(i =>
		//    {
		//        double m, v;
		//        i.GetMeanAndVariance(out m, out v);
		//        return m + .5 * v;
		//    }
		//        ));
		//}
	}


	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	public static class SoftmaxOp_Bohning
	{
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			// this formula is based on the following bound, which resembles a second-order Taylor approximation:
			// E[log(sum_k exp(x_k))] <= log(sum_k exp(m_k)) + 0.5 trace(V H)
			// where V is the (diagonal) covariance matrix of x and H is Bohning's bound on the Hessian matrix of log-sum-exp
			// H = 0.5*(I - 11'/K)
			double lse = double.NegativeInfinity;
			double trace = 0;
			double avgLog = 0;
			int K = x.Count;
			for (int i = 0; i < K; i++) {
				double m, v;
				x[i].GetMeanAndVariance(out m, out v);
				lse = MMath.LogSumExp(lse, m);
				trace += v;
				avgLog += (softmax.PseudoCount[i] - 1)*m;
			}
			trace *= 0.5*(1 - 1.0/K);
			avgLog -= (softmax.TotalCount - K)*(0.5*trace + lse);
			return avgLog -softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax);
		}

		//public static double ExpectationLogSumExp<GaussianList>(GaussianList x)
		//    where GaussianList : IList<Gaussian>
		//{
		//    double lse = double.NegativeInfinity;
		//    double trace = 0;
		//    int K = x.Count;
		//    for (int i = 0; i < K; i++) {
		//        double m, v;
		//        x[i].GetMeanAndVariance(out m, out v);
		//        lse = MMath.LogSumExp(lse, m);
		//        trace += v;
		//    }
		//    trace *= 0.5*(1 - 1.0/K);
		//    return .5*trace + lse;
		//}

		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			double lse = double.NegativeInfinity;
			double trace = 0;
			int K = x.Count;
			Vector meanLog = Vector.Zero(K);
			for (int i = 0; i < K; i++) {
				double m, v;
				x[i].GetMeanAndVariance(out m, out v);
				lse = MMath.LogSumExp(lse, m);
				trace += v;
				meanLog[i] = m;
			}
			trace *= 0.5*(1 - 1.0/K);
			double denom = 0.5*trace+lse;
			for (int i = 0; i < K; i++) {
				meanLog[i] -= denom;
			}
			result.SetMeanLog(meanLog);
			return result;
		}
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [SkipIfAnyUniform] IList<Gaussian> x, GaussianList result)
						where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			double c = softmax.TotalCount - K;
			double prec = c*0.5*(1 - 1.0/K);
			double lse = double.NegativeInfinity;
			for (int i = 0; i < K; i++) {
				double m = x[i].GetMean();
				lse = MMath.LogSumExp(lse, m);
			}
			for (int i = 0; i < K; i++) {
				double m = x[i].GetMean();
				double s = Math.Exp(m - lse);
				double meanTimesPrec = m*prec - c*s + (softmax.PseudoCount[i]-1);
				result[i] = Gaussian.FromNatural(meanTimesPrec, prec);
			}
			return result;
		}
	}

	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	public static class SoftmaxOp_Taylor
	{
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			// this formula is based on the following second-order Taylor approximation:
			// E[log(sum_k exp(x_k))] =approx log(sum_k exp(m_k)) + 0.5 trace(V H)
			// where V is the (diagonal) covariance matrix of x and H is the Hessian matrix of log-sum-exp
			double lse = double.NegativeInfinity;
			double trace = 0;
			double avgLog = 0;
			int K = x.Count;
			for (int i = 0; i < K; i++) {
				double m = x[i].GetMean();
				lse = MMath.LogSumExp(lse, m);
			}
			for (int i = 0; i < K; i++) {
				double m, v;
				x[i].GetMeanAndVariance(out m, out v);
				double s = Math.Exp(m - lse);
				trace += v*s*(1-s);
				avgLog += (softmax.PseudoCount[i] - 1)*m;
			}
			avgLog -= (softmax.TotalCount - K)*(0.5*trace + lse);
			return avgLog -softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax);
		}

		//public static double ExpectationLogSumExp<GaussianList>(GaussianList x)
		//        where GaussianList : IList<Gaussian>
		//{
		//    double lse = double.NegativeInfinity;
		//    double trace = 0;
		//    int K = x.Count;
		//    for (int i = 0; i < K; i++)
		//    {
		//        double m = x[i].GetMean();
		//        lse = MMath.LogSumExp(lse, m);
		//    }
		//    for (int i = 0; i < K; i++)
		//    {
		//        double m, v;
		//        x[i].GetMeanAndVariance(out m, out v);
		//        double s = Math.Exp(m - lse);
		//        trace += .5 * v * s * (1 - s);
		//    }
		//    return trace + lse;
		//}


		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			double lse = double.NegativeInfinity;
			double trace = 0;
			int K = x.Count;
			Vector meanLog = Vector.Zero(K);
			for (int i = 0; i < K; i++) {
				double m = x[i].GetMean();
				lse = MMath.LogSumExp(lse, m);
			}
			for (int i = 0; i < K; i++) {
				double m, v;
				x[i].GetMeanAndVariance(out m, out v);
				double s = Math.Exp(m - lse);
				trace += v*s*(1-s);
				meanLog[i] = m;
			}
			double denom = 0.5*trace+lse;
			for (int i = 0; i < K; i++) {
				meanLog[i] -= denom;
			}
			result.SetMeanLog(meanLog);
			return result;
		}
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [SkipIfAnyUniform] IList<Gaussian> x, GaussianList to_x)
						where GaussianList : IList<Gaussian>
		{
			GaussianList result = to_x;
			int K = x.Count;
			double c = softmax.TotalCount - K;
			double lse = double.NegativeInfinity;
			for (int i = 0; i < K; i++) {
				double m = x[i].GetMean();
				lse = MMath.LogSumExp(lse, m);
			}
			double step = 0.5;
			for (int k = 0; k < K; k++) {
				double m = x[k].GetMean();
				double s = Math.Exp(m - lse);
				double prec = c*s*(1-s);
				double meanTimesPrec = m*prec - c*s + (softmax.PseudoCount[k]-1);
				Gaussian rk = Gaussian.FromNatural(meanTimesPrec, prec);
				rk.Precision = step * rk.Precision + (1 - step) * to_x[k].Precision;
				rk.MeanTimesPrecision = step * rk.MeanTimesPrecision + (1 - step) * to_x[k].MeanTimesPrecision;
				result[k] = rk;
			}
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the bound from Knowles and Minka (2011), followed by 
	/// nonconjugate VMP. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>), Default=true)]
	[Quality(QualityBand.Preview)]
	[Buffers("A", "history")]
	public static class SoftmaxOp_KM11
	{
		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static Vector AInit([IgnoreDependency] IList<Gaussian> x)
		{
			return Vector.Constant(x.Count, 1.0 / x.Count);
		}

		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Vector A([Proper] IList<Gaussian> x, Vector a)
		{
			return SoftmaxOp_KM11_LBFGS.A(x, a);
		}

		public static List<IList<Gaussian>> HistoryInit()
		{
			return new List<IList<Gaussian>>();
		}
		public static List<IList<Gaussian>> History([Proper] IList<Gaussian> x, List<IList<Gaussian>> history)
		{
			history.Add((GaussianArray)((GaussianArray)x).Clone());
			return history;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, [SkipIfUniform] Vector a, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			return SoftmaxOp_KM11_LBFGS.SoftmaxAverageLogarithm(x, a, result);
		}
		[Skip]
		public static Dirichlet SoftmaxAverageLogarithmInit([IgnoreDependency] IList<Gaussian> x)
		{
			return Dirichlet.Uniform(x.Count);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Vector a, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			return SoftmaxOp_KM11_LBFGS.AverageLogFactor(x, a, softmax, to_softmax);
		}

		public static bool useBounds = false;

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_x">Previous outgoing message to 'x'.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>The outgoing VMP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [SkipIfAnyUniform] IList<Gaussian> x, GaussianList to_x, Vector a)
						where GaussianList : IList<Gaussian>
		{
			GaussianList result = to_x;
			int K = x.Count;
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var counts = softmax.PseudoCount - 1;
			double sum_n = softmax.TotalCount - K;

			double step = (LogisticOp_SJ99.global_step == 0.0) ? 1.0 : (Rand.Double() * LogisticOp_SJ99.global_step); // random damping helps convergence, especially with parallel updates
			//double step = 1 ;
			//double step = Rand.Double();
			if (false) {
				// assumes "a" is optimized
				for (int k = 0; k < K; k++) {
					var rk = result[k];
					rk.Precision = sum_n * a[k] * (1.0 - a[k]);
					rk.MeanTimesPrecision = m[k] * rk.Precision + counts[k] - sum_n * a[k];
					if (step != 1.0) {
						rk.Precision = step * rk.Precision + (1 - step) * to_x[k].Precision;
						rk.MeanTimesPrecision = step * rk.MeanTimesPrecision + (1 - step) * to_x[k].MeanTimesPrecision;
					}
					result[k] = rk;
				}
			} else if (false) {
				// use history to get better updates
				double[] dm, dv;
				GetDerivatives(counts, sum_n, x, a, out dm, out dv);
				if (true) {
					double[] dm2, dv2;
					//IList<Gaussian> x2 = history[history.Count-2];
					int j = Rand.Int(K);
					//Console.WriteLine(j);
					//if (j == 16) Console.WriteLine(j);
					for (int k = 0; k < K; k++) {
						if (k != j) continue;
						var rk = result[k];
						double m2;// = x2[k].GetMean();
						m2 = (m[k] * dv[k] + dm[k]) / dv[k];
						GaussianArray xt = (GaussianArray)((GaussianArray)x).Clone();
						xt[k] = Gaussian.FromMeanAndPrecision(m2, dv[k]);
						GetDerivatives(counts, sum_n, xt, a, out dm2, out dv2);
						rk.Precision = step * dv[k] + (1 - step) * dv2[k];
						rk.Precision = dv[k];
						// estimate 2nd derivative using finite differences
						//rk.Precision = -(dm[k] - dm2[k])/(m[k] - m2);
						if (rk.Precision <= 0 || Double.IsNaN(rk.Precision)) throw new Exception("stop");
						rk.MeanTimesPrecision = step * (m[k] * rk.Precision + dm[k]) + (1 - step) * (m2 * rk.Precision + dm2[k]);
						result[k] = rk;
					}
				} else {
					for (int k = 0; k < K; k++) {
						var rk = result[k];
						rk.Precision = dv[k];
						rk.MeanTimesPrecision = m[k] * rk.Precision + dm[k];
						result[k] = rk;
					}
				}
			} else {
				// does not assume "a" is optimized
				double logSumExp = SoftmaxOp_KM11_LBFGS.ExpectationLogSumExp_Helper(Vector.FromArray(m), Vector.FromArray(v), a);

				double precBound = 3;
				double meanBound = 3;
				bool showBounds = false;
				for (int k = 0; k < K; k++) {
					double logistic_k = Math.Exp(m[k] + (1.0 - 2.0 * a[k]) * v[k] / 2.0 - logSumExp); //= a[k]
					var rk = result[k];
					rk.Precision = sum_n * ((1.0 - 2.0 * a[k]) * logistic_k + a[k] * a[k]);
					if (useBounds && false) {
						double oldPrec = to_x[k].Precision;
						double xPrec = x[k].Precision;
						double xPrecNew = xPrec + (rk.Precision - oldPrec);
						double maxPrec = xPrec * precBound;
						double minPrec = xPrec / precBound;
						if (xPrecNew > maxPrec) {
							if (showBounds) Console.WriteLine("maxPrec k={0}", k);
							// set xPrecNew = maxPrec
							rk.Precision = maxPrec - xPrec + oldPrec;
						} else if (xPrecNew < minPrec) {
							if (showBounds) Console.WriteLine("minPrec k={0}", k);
							// set xPrecNew = minPrec
							rk.Precision = minPrec - xPrec + oldPrec;
						}
						xPrecNew = xPrec + (rk.Precision - oldPrec);
						//double precDiff = Math.Exp(Math.Abs(Math.Log(xPrecNew / xPrec)));
					}
					if (false && !to_x[k].IsUniform()) {
						// damp Precision before computing MeanTimesPrecision
						// doesn't help as much as damping at end
						rk.Precision = step * rk.Precision + (1 - step) * to_x[k].Precision;
					}
					rk.MeanTimesPrecision = m[k] * rk.Precision + counts[k] - sum_n * logistic_k;
					if (useBounds && true) {
						double oldPrec = to_x[k].Precision;
						double xPrec = x[k].Precision;
						double xPrecNew = xPrec + (rk.Precision - oldPrec);
						double oldMeanTimesPrecision = to_x[k].MeanTimesPrecision;
						double maxMean = m[k] + meanBound / Math.Sqrt(xPrec);
						double minMean = m[k] - meanBound / Math.Sqrt(xPrec);
						double xMeanNew = (x[k].MeanTimesPrecision + rk.MeanTimesPrecision - oldMeanTimesPrecision) / xPrecNew;
						if (xMeanNew > maxMean) {
							xMeanNew = maxMean;
							if (showBounds) Console.WriteLine("maxMean k={0}", k);
						} else if (xMeanNew < minMean) {
							xMeanNew = minMean;
							if (showBounds) Console.WriteLine("minMean k={0}", k);
						}
						rk.MeanTimesPrecision = xMeanNew * xPrecNew - x[k].MeanTimesPrecision + oldMeanTimesPrecision;
						//double meanDiff = Math.Abs(xMeanNew - m[k])*Math.Sqrt(xPrec);
					}
					if (true && step != 1.0) {
						rk.Precision = step * rk.Precision + (1 - step) * to_x[k].Precision;
						rk.MeanTimesPrecision = step * rk.MeanTimesPrecision + (1 - step) * to_x[k].MeanTimesPrecision;
					}
					if (false) {
						// this update should help, but has little effect
						// it seems that interaction between the x's is not the main issue
						Gaussian xkNew = x[k] / to_x[k] * rk;
						xkNew.GetMeanAndVariance(out m[k], out v[k]);
						GaussianArray xNew = (GaussianArray)((GaussianArray)x).Clone();
						xNew[k] = xkNew;
						x = xNew;
						a = A(x, a);
						logSumExp = SoftmaxOp_KM11_LBFGS.ExpectationLogSumExp_Helper(Vector.FromArray(m), Vector.FromArray(v), a);
					}
					result[k] = rk;
				}
			}
			return result;
		}

		private static void GetDerivatives(Vector counts, double sum_n, IList<Gaussian> x, Vector a, out double[] dm, out double[] dv)
		{
			int K = x.Count;
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var ms = Vector.FromArray(m);
			var vs = Vector.FromArray(v);
			dm = new double[K];
			dv = new double[K];
			double logSumExp = SoftmaxOp_KM11_LBFGS.ExpectationLogSumExp_Helper(ms, vs, a);
			for (int k = 0; k < K; k++) {
				double logistic_k = Math.Exp(m[k] + (1.0 - 2.0 * a[k]) * v[k] / 2.0 - logSumExp); //= a[k]
				dv[k] = sum_n * ((1.0 - 2.0 * a[k]) * logistic_k + a[k] * a[k]);
				dm[k] = counts[k] - sum_n * logistic_k;
			}
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses a generalisation of the tilted bound used in Saul Jordan 1999, followed by 
	/// nonconjugate VMP. This approach is linear in the dimension K. 
	/// 
	/// This will replace SaulJordanSoftmaxOp_NCVMP once the functionality is confirmed to be equivalent. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	[Buffers("Asj", "Abouchard")]
	public static class AdaptiveSoftmaxOp_NCVMP_Sparse
	{
		public static IList<double> AsjInit([IgnoreDependency] IList<Gaussian> x)
		{
			//return Vector.Constant(x.Count, 1.0 / x.Count, x.IsSparse() ? Sparsity.Sparse : Sparsity.Dense);
			return x.ListSelect(i => 1.0 / x.Count);
		}

		public static IList<double> Asj([Proper] IList<Gaussian> x, IList<double> asj)
		{
			return SoftmaxOp_KM11_LBFGS_Sparse.A(x, asj);
		}

		public static double AbouchardInit()
		{
			return .5;
		}

		public static double Abouchard<GaussianList>([SkipIfAnyUniform] GaussianList x, double abouchard) where GaussianList : IList<Gaussian>
		{
			return BouchardSoftmaxOp_Sparse.A(x, abouchard);
		}


		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, [SkipIfUniform] IList<double> asj, double abouchard, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			if (SoftmaxOp_KM11_Sparse.AverageLogFactor(x, asj, softmax, to_softmax)
                > BouchardSoftmaxOp_Sparse.AverageLogFactor(x, softmax, to_softmax, abouchard))
				return SoftmaxOp_KM11_Sparse.SoftmaxAverageLogarithm(x, asj, to_softmax);
			else
				return BouchardSoftmaxOp_Sparse.SoftmaxAverageLogarithm(x, abouchard, to_softmax);
		}

		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, [SkipIfUniform] IList<double> asj, double abouchard, Dirichlet softmax, Dirichlet to_softmax, GaussianList to_x)
				where GaussianList : IList<Gaussian>
		{
			if (SoftmaxOp_KM11_Sparse.AverageLogFactor(x, asj, softmax, to_softmax)
                > BouchardSoftmaxOp_Sparse.AverageLogFactor(x, softmax, to_softmax, abouchard))
				return SoftmaxOp_KM11_Sparse.XAverageLogarithm(softmax, x, to_x, asj);
			else
				return BouchardSoftmaxOp_Sparse.XAverageLogarithm(softmax, x, abouchard, to_x);
		}

		public static double LogAverageFactor<GaussianList>([SkipIfAllUniform] GaussianList x, [SkipIfUniform] IList<double> asj, double abouchard, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			double SJbound = SoftmaxOp_KM11_Sparse.AverageLogFactor(x, asj, softmax, to_softmax);
			double BouchardBound = BouchardSoftmaxOp_Sparse.AverageLogFactor(x, softmax, to_softmax, abouchard);
			return Math.Max(SJbound, BouchardBound);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the bound in Knowles and Minka (2011), followed by 
	/// nonconjugate VMP. This approach is linear in the dimension K. 
	/// 
	/// This will replace SaulJordanSoftmaxOp_NCVMP once the functionality is confirmed to be equivalent. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	[Buffers("A")]
	public static class SoftmaxOp_KM11_Sparse
	{

		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static IList<double> AInit([IgnoreDependency] IList<Gaussian> x)
		{
			//return Vector.Constant(x.Count, 1.0 / x.Count, x.IsSparse() ? Sparsity.Sparse : Sparsity.Dense);
			return x.ListSelect(i => 1.0 / x.Count);
		}

		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static IList<double> A([Proper] IList<Gaussian> x, IList<double> a)
		{
			Func<Gaussian, double, double> f = (xi, ai) => .5 - (Math.Log(ai) - xi.GetMean()) / xi.GetVariance();
			var a1 = f.Map(x, a);
			double C = (a1.Sum() - 1.0) / x.Select(o => o.Precision).Sum();
			Func<Gaussian, double, double> g = (xi, ai) => ai - C * xi.Precision;
			IList<double> newa = a.IsSparse() ? (IList<double>)SparseList<double>.FromSize(a.Count) : new List<double>(a.Count);
			newa.SetTo(g.Map(x, a1));

			var normala = MMath.Softmax(x.ListZip(a, (i, j) => {
				double m, v;
				i.GetMeanAndVariance(out m, out v);
				double res = m + (1.0 - 2.0 * j) * v / 2.0;
				if (double.IsNaN(res))
					throw new ApplicationException("nans");
				return res;
			}));

			if (true) {
				double step = 0.99;
				Func<double,double, double> damp = (oldv, newv) => (oldv*(1-step) + newv*step);
				a.SetTo(damp.Map(a, normala));
				return a;
			}

			if (logSumExpFull(x, normala) > logSumExpFull(x, newa))
				return newa;
			else
				return normala;
		}


		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, [SkipIfUniform] IList<double> a, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			Func<Gaussian, double, double> f = (xi, ai) => xi.GetVariance() * ai * ai;
			double sum = .5 * f.Map(x, a).EnumerableSum(i => i);
			f = (xi, ai) =>
					xi.GetMean() + (1.0 - 2.0 * ai) * 0.5 * xi.GetVariance();
			sum += MMath.LogSumExpSparse(f.Map(x, a));
			result.SetMeanLog(Vector.FromList(x.ListSelect(xi => xi.GetMean() - sum)));
			if (result.PseudoCount.Any(o => double.IsNaN(o)))
				throw new ApplicationException("nans");
			return result;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, IList<double> a, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			double sum_n = softmax.PseudoCount.EnumerableSum(i => i - 1.0);

			Func<Gaussian, double, double> f = (i, j) => i.GetMean() * (j - 1.0);
			double innerProduct = f.Map(x, softmax.PseudoCount).EnumerableSum(i => i);


			return innerProduct - sum_n * logSumExpFull(x, a) - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax);
		}

		public static double logSumExpFull(IEnumerable<Gaussian> x, IEnumerable<double> a)
		{
			double logSumExp = ExpectationLogSumExp_Helper(x, a);
			Func<Gaussian, double, double> f = (xi, ai) => xi.GetVariance() * ai * ai;
			double halfSumA2v = .5 * f.Map(x, a).EnumerableSum(i => i);
			return logSumExp + halfSumA2v;
		}


		internal static double ExpectationLogSumExp_Helper(IEnumerable<Gaussian> x, IEnumerable<double> a)
		{
			// returns logsumexp(m + v*(0.5 - a))
			Func<Gaussian, double, double> f = (xi, ai) => {
				double m, v;
				xi.GetMeanAndVariance(out m, out v);
				return m + v * (.5 - ai);
			};
			return MMath.LogSumExp(f.Map(x, a));
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_x">Previous outgoing message to 'x'.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>The outgoing VMP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [SkipIfAnyUniform] GaussianList x, GaussianList to_x, IList<double> a)
				where GaussianList : IList<Gaussian>
		{
			double sum_n = softmax.PseudoCount.EnumerableSum(i => i - 1.0);
			double step = (LogisticOp_SJ99.global_step == 0.0) ? 1.0 : (LogisticOp_SJ99.global_step * Rand.Double());  // random damping helps convergence, especially with parallel updates
			// does not assume "a" is optimized
			double logSumExp = ExpectationLogSumExp_Helper(x, a);
			Func<Gaussian, double, double, Gaussian, Gaussian> f = (xk, ak, countsk, to_xk) => {
				double mk, vk;
				xk.GetMeanAndVariance(out mk, out vk);
				double logistic_k = Math.Exp(mk + (1.0 - 2.0 * ak) * vk / 2.0 - logSumExp); //= ak
				var rk = new Gaussian(); // ideally would get this from result (i.e. to_x)
				rk.Precision = sum_n * ((1.0 - 2.0 * ak) * logistic_k + ak * ak);
				rk.MeanTimesPrecision = mk * rk.Precision + countsk - 1.0 - sum_n * logistic_k;
				if (step != 1.0) {
					rk.Precision = step * rk.Precision + (1 - step) * to_xk.Precision;
					rk.MeanTimesPrecision = step * rk.MeanTimesPrecision + (1 - step) * to_xk.MeanTimesPrecision;
				}
				if (double.IsNaN(rk.Precision) || double.IsNaN(rk.MeanTimesPrecision))
					throw new ApplicationException("nans");
				return rk;
			};
			to_x.SetTo(f.Map(x, a, softmax.PseudoCount, to_x));
			return to_x;
		}
	}


	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the simple first order Taylor series expansion from Blei and Lafferty (2006), followed by
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	public static class SoftmaxOp_BL06
	{

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			Func<Gaussian, double> f = (xi) =>
					xi.GetMean() + 0.5 * xi.GetVariance();
			double sum = MMath.LogSumExpSparse(f.Map(x));
			result.SetMeanLog(Vector.FromList(x.ListSelect(xi => xi.GetMean() - sum)));
			return result;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			double sum_n = softmax.PseudoCount.EnumerableSum(i => i - 1.0);
			double logSumExp = ExpectationLogSumExp_Helper(x);
			Func<Gaussian, double, double> f = (i, j) => i.GetMean() * (j - 1.0);
			double innerProduct = f.Map(x, softmax.PseudoCount).EnumerableSum(i => i);
			return innerProduct - sum_n * logSumExp - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax);
		}


		public static double ExpectationLogSumExp_Helper(IList<Gaussian> x)
		{
			// returns logsumexp(m + v*0.5)
			Func<Gaussian, double> f = (xi) => {
				double m, v;
				xi.GetMeanAndVariance(out m, out v);
				return m + v * .5;
			};
			return MMath.LogSumExp(f.Map(x));
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_x">Previous outgoing message to 'x'.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>The outgoing VMP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [SkipIfUniform] GaussianList x, GaussianList to_x)
				where GaussianList : IList<Gaussian>
		{
			double sum_n = softmax.PseudoCount.EnumerableSum(i => i - 1.0);
			double step = Rand.Double() * 0.5; // random damping helps convergence, especially with parallel updates
			double logSumExp = ExpectationLogSumExp_Helper(x);
			Func<Gaussian, double, Gaussian, Gaussian> f = (xk, countsk, to_xk) => {
				double mk, vk;
				xk.GetMeanAndVariance(out mk, out vk);
				double logistic_k = Math.Exp(mk + vk / 2.0 - logSumExp); //= ak
				var rk = new Gaussian(); // ideally would get this from result (i.e. to_x)
				rk.Precision = sum_n * logistic_k;
				rk.MeanTimesPrecision = mk * rk.Precision + countsk - 1.0 - sum_n * logistic_k;
				if (step != 1.0) {
					rk.Precision = step * rk.Precision + (1 - step) * to_xk.Precision;
					rk.MeanTimesPrecision = step * rk.MeanTimesPrecision + (1 - step) * to_xk.MeanTimesPrecision;
				}
				return rk;
			};
			to_x.SetTo(f.Map(x, softmax.PseudoCount, to_x));
			return to_x;
		}
	}


	/// <summary>
	/// Provides outgoing messages for the vector version of <see cref="MMath.Softmax(Vector)"/>, given random arguments to the function.
	/// This implementation uses the bound in Knowles and Minka (2011), followed by 
	/// nonconjugate VMP. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(Vector))]
	[Quality(QualityBand.Preview)]
	[Buffers("A")]
	public static class VectorSoftmaxOp_KM11
	{

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(softmax,x))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static Vector AInit([IgnoreDependency] VectorGaussian x)
		{
			return Vector.Constant(x.Dimension, 1.0 / x.Dimension);
		}

		public static Gaussian[] VectorGaussianToGaussianList(VectorGaussian x)
		{
			var x2 = new Gaussian[x.Dimension];
			var m = Vector.Zero(x.Dimension);
			var v = new PositiveDefiniteMatrix(x.Dimension, x.Dimension);
			x.GetMeanAndVariance(m, v);
			for (int i = 0; i < x.Dimension; i++) {
				x2[i] = Gaussian.FromMeanAndVariance(m[i], v[i, i]);
			}
			return x2;
		}


		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Vector A([Proper, SkipIfUniform] VectorGaussian x, Vector a)
		{

			return SoftmaxOp_KM11_LBFGS.A(VectorGaussianToGaussianList(x), a);
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm([SkipIfUniform] VectorGaussian x, [SkipIfUniform] Vector a, Dirichlet result)
		{
			return SoftmaxOp_KM11_LBFGS.SoftmaxAverageLogarithm(VectorGaussianToGaussianList(x), a, result);
		}
		[Skip]
		public static Dirichlet SoftmaxAverageLogarithmInit([IgnoreDependency] VectorGaussian x)
		{
			return Dirichlet.Uniform(x.Dimension);
		}


		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_x">Previous outgoing message to 'x'.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>The outgoing VMP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static VectorGaussian XAverageLogarithm([SkipIfUniform] Dirichlet softmax, [SkipIfUniform] VectorGaussian x, VectorGaussian to_x, Vector a)
		{
			VectorGaussian result = to_x;
			int K = x.Dimension;
			var m = Vector.Zero(x.Dimension);
			var v = new PositiveDefiniteMatrix(x.Dimension, x.Dimension);
			x.GetMeanAndVariance(m, v);
			var counts = softmax.PseudoCount - 1;
			double sum_n = counts.Sum();
			for (int k = 0; k < K; k++) {
				result.Precision[k, k] = sum_n * a[k] * (1.0 - a[k]);
				result.MeanTimesPrecision[k] = counts[k] - sum_n * a[k];
				for (int l = 0; l < K; l++)
					if (l != k) {
						result.Precision[k, l] = -sum_n * a[k] * a[l] * (a[l] + a[k]);
					}
			}
			result.MeanTimesPrecision.SetToSum(result.MeanTimesPrecision, result.Precision * m);
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the bound in Knowles and Minka (2011), followed by 
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	[Buffers("A")]
	public static class SoftmaxOp_KM11_LBFGS
	{

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, Vector a, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			Vector ms, vs;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out ms, out vs);
			var ns = softmax.PseudoCount - 1.0;
			double sum_n = ns.Sum();
			double logSumExp = ExpectationLogSumExp_Helper(ms, vs, a);
			double halfSumA2v = .5 * (a * a * vs).Sum();
			return ns.Inner(ms) - sum_n * (logSumExp + halfSumA2v) - softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax);
		}

		internal static double ExpectationLogSumExp_Helper(Vector ms, Vector vs, Vector a)
		{
			// returns logsumexp(m + v*(0.5 - a))
			var temp = Vector.Zero(a.Count, ms.Sparsity);
			temp.SetToProduct(a, -1.0);
			temp.SetToSum(temp, .5);
			temp.SetToProduct(temp, vs);
			temp.SetToSum(temp, ms);
			return MMath.LogSumExp(temp);
		}

		/// <summary>
		/// Function to evaluate this factor's KL divergence contribution, and gradient (if grad is not null). 
		/// </summary>
		/// <param name="mu">Prior means</param>
		/// <param name="s2">Prior variances</param>
		/// <param name="a">Variational parameter vector a</param>
		/// <param name="x">x[1..K]: posterior mean, x[K+1..2K]: posterior log variance</param>
		/// <param name="ns">Dirichlet counts-1</param>
		/// <param name="grad">Vector to store the gradient in</param>
		/// <returns>KL divergence</returns>
		public static double GradientAndValueAtPoint(double[] mu, double[] s2, Vector a, Vector x, Vector ns, Vector grad)
		{
			int K = x.Count / 2;
			var ms = Vector.Zero(K);
			var vs = Vector.Zero(K);
			for (int k = 0; k < K; k++) {
				ms[k] = x[k];
				vs[k] = Math.Exp(x[k + K]);
			}
			double sum_n = ns.Sum();
			double logSumExp = ExpectationLogSumExp_Helper(ms, vs, a);
			double halfSumA2v = .5 * (a * a * vs).Sum();
			double kl_value = sum_n * (logSumExp + halfSumA2v);
			for (int k = 0; k < K; k++) {
				if (s2[k] > 0.0)
					kl_value += -.5 * (1.0 + x[K + k]) + .5 * (vs[k] + (ms[k] - mu[k]) * (ms[k] - mu[k])) / s2[k] + .5 * Math.Log(s2[k]) - ns[k] * ms[k];
				if (double.IsInfinity(kl_value) || double.IsNaN(kl_value))
					throw new ApplicationException("KL became undefined while optimising variational parameters of the softmax factor. Try using Blei06SoftmaxOp.");
				if (grad != null) {
					if (s2[k] == 0.0) {
						grad[k] = 0.0;
						grad[K + k] = 0.0;
					} else {
						double logistic_k = Math.Exp(ms[k] + (1.0 - 2.0 * a[k]) * vs[k] / 2.0 - logSumExp);
						grad[k] = (ms[k] - mu[k]) / s2[k] - ns[k] + sum_n * logistic_k;
						grad[K + k] = -.5 + .5 * vs[k] / s2[k] + .5 * vs[k] * sum_n * ((1.0 - 2.0 * a[k]) * logistic_k + a[k] * a[k]);
					}
				}
			}

			return kl_value;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, Vector a, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			Vector meanLog = Vector.Zero(K);
			double[] exponent = new double[K];
			double sum = 0.0;
			for (int i = 0; i < K; i++) {
				double m, v;
				x[i].GetMeanAndVariance(out m, out v);
				meanLog[i] = m;
				exponent[i] = m + (1.0 - 2.0 * a[i]) * 0.5 * v;
				sum += v * a[i] * a[i];
			}
			sum = 0.5 * sum + MMath.LogSumExp(exponent);
			for (int i = 0; i < K; i++) {
				meanLog[i] -= sum;
			}
			result.SetMeanLog(meanLog);
			return result;
		}
		[Skip]
		public static Dirichlet SoftmaxAverageLogarithmInit([IgnoreDependency] IList<Gaussian> x)
		{
			return Dirichlet.Uniform(x.Count);
		}

		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static Vector AInit([IgnoreDependency] IList<Gaussian> x)
		{
			return Vector.Constant(x.Count, 1.0 / x.Count, x.IsSparse() ? Sparsity.Sparse : Sparsity.Dense);
		}

		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Vector A([Proper] IList<Gaussian> x, Vector a)
		{
			Vector m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var temp = Vector.Zero(m.Count, m.Sparsity);
			// temp = m + (1-2a)*v/2
			temp.SetToProduct(a, -1.0);
			temp.SetToSum(temp, .5);
			temp.SetToProduct(temp, v);
			temp.SetToSum(temp, m);
			a = MMath.Softmax(temp);
			return a;
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [Stochastic, SkipIfAllUniform] IList<Gaussian> x, Vector a, GaussianList result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			var prior = new Gaussian[K];
			for (int k = 0; k < K; k++) {
				prior[k] = x[k] / result[k];
			}
			double[] mu, s2;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(prior, out mu, out s2);
			double[] m, v;
			SoftmaxOp_BL06_LBFGS.GetMeanAndVariance(x, out m, out v);
			var counts = softmax.PseudoCount - 1;
			int evalCounter = 0;
			// z[1..K]: posterior means, z[K+1..2K]: log posterior variances
			var z = Vector.Zero(K * 2);
			for (int k = 0; k < K; k++) {
				z[k] = m[k];
				z[K + k] = Math.Log(v[k]);
			}
			double start0 = GradientAndValueAtPoint(mu, s2, Vector.Zero(K), z, counts, null);
			double startingValue = GradientAndValueAtPoint(mu, s2, a, z, counts, null);
			var s = new LBFGS(5);
			s.MaximumStep = 1e3;
			s.MaximumIterations = 100;
			s.Epsilon = 1e-10;
			s.convergenceCriteria = BFGS.ConvergenceCriteria.Objective;
			z = s.Run(z, 1.0, delegate(Vector y, ref Vector grad) { evalCounter++; return GradientAndValueAtPoint(mu, s2, a, y, counts, grad); });
			for (int k = 0; k < K; k++) {
				m[k] = z[k];
				v[k] = Math.Exp(z[K + k]);
				var rk = result[k];
				if (!prior[k].IsPointMass)
					rk.SetToRatio(Gaussian.FromMeanAndVariance(m[k], v[k]), prior[k]);
				result[k] = rk;
			}
			double endValue = GradientAndValueAtPoint(mu, s2, a, z, counts, null);
			//Console.WriteLine("Went from {0} to {1} in {2} steps, {3} evals", startingValue, endValue, s.IterationsPerformed, evalCounter);
			if (startingValue < endValue)
				Console.WriteLine("Warning: LBFGS resulted in an increased objective function");
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the bound in Knowles and Minka (2011), followed by 
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	[Buffers("A")]
	public static class SoftmaxOp_KM11_LBFGS_Sparse
	{

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		public static double AverageLogFactor<GaussianList>([SkipIfAnyUniform] GaussianList x, IList<double> a, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			return SoftmaxOp_KM11_Sparse.AverageLogFactor(x, a, softmax, to_softmax);
		}

		//internal static double ExpectationLogSumExp_Helper(Vector ms, Vector vs, Vector a)
		//{
		//    // returns logsumexp(m + v*(0.5 - a))
		//    var temp = Vector.Zero(a.Count, ms.Sparsity);
		//    temp.SetToProduct(a, -1.0);
		//    temp.SetToSum(temp, .5);
		//    temp.SetToProduct(temp, vs);
		//    temp.SetToSum(temp, ms);
		//    return MMath.LogSumExp(temp);
		//}

		/// <summary>
		/// Function to evaluate this factor's KL divergence contribution, and gradient (if grad is not null). 
		/// </summary>
		/// <param name="mu">Prior means</param>
		/// <param name="s2">Prior variances</param>
		/// <param name="a">Variational parameter vector a</param>
		/// <param name="x">x[1..K]: posterior mean, x[K+1..2K]: posterior log variance</param>
		/// <param name="ns">Dirichlet counts-1</param>
		/// <param name="grad">Vector to store the gradient in</param>
		/// <returns>KL divergence</returns>
		public static double GradientAndValueAtPoint(IList<Gaussian> prior, IList<double> a, Vector xm, Vector lxv, Vector ns, Vector[] grad)
		{
			// remember to subtract one from ns !
			int K = xm.Count;
			double sum_n = ns.Sum(o => o - 1);
			double logSumExp = MMath.LogSumExpSparse(new Func<double, double, double, double>((m, lv, ak) => m + Math.Exp(lv) * (.5 - ak)).Map(xm, lxv, a));
			double halfSumA2v = .5 * new Func<double, double, double>((ak, lv) => ak * ak * Math.Exp(lv)).Map(a, lxv).Sum();

			if (grad != null) {
				Func<Gaussian, double, double, double, double, double> gm = (priork, msk, lvsk, ak, nsk) => {
					double muk, s2k;
					priork.GetMeanAndVariance(out muk, out s2k);
					if (s2k == 0.0)
						return 0.0;
					else {
						double logistic_k = Math.Exp(msk + (1.0 - 2.0 * ak) * Math.Exp(lvsk) / 2.0 - logSumExp);
						return (msk - muk) / s2k - (nsk - 1) + sum_n * logistic_k;
					}
				};

				grad[0].SetTo(gm.Map(prior, xm, lxv, a, ns));

				Func<Gaussian, double, double, double, double, double> gl = (priork, msk, lvsk, ak, nsk) =>
				{
					double muk, s2k;
					priork.GetMeanAndVariance(out muk, out s2k);
					if (s2k == 0.0)
						return 0.0;
					else {
						double logistic_k = Math.Exp(msk + (1.0 - 2.0 * ak) * Math.Exp(lvsk) / 2.0 - logSumExp);
						double vsk = Math.Exp(lvsk);
						return -.5 + .5 * vsk / s2k + .5 * vsk * sum_n * ((1.0 - 2.0 * ak) * logistic_k + ak * ak);
					}
				};

				grad[1].SetTo(gl.Map(prior, xm, lxv, a, ns));
			}

			Func<Gaussian, double, double, double, double> f = (priork, msk, lxvk, nsk) => {
				double muk, s2k;
				priork.GetMeanAndVariance(out muk, out s2k);
				return s2k > 0.0 ?
                     -.5 * (1.0 + lxvk) + .5 * (Math.Exp(lxvk) + (msk - muk) * (msk - muk)) / s2k + .5 * Math.Log(s2k) - (nsk - 1) * msk
                     : 0.0;
			};

			var res = f.Map(prior, xm, lxv, ns).Sum() + sum_n * (logSumExp + halfSumA2v);
			return res;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, IList<double> a, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			return SoftmaxOp_KM11_Sparse.SoftmaxAverageLogarithm(x, a, result);
		}

		/// <summary>
		/// Initialise the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Initial value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static IList<double> AInit([IgnoreDependency] IList<Gaussian> x)
		{
			return SoftmaxOp_KM11_Sparse.AInit(x);
		}

		/// <summary>
		/// Update the buffer 'A'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <returns>New value of buffer 'A'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static IList<double> A([Proper] IList<Gaussian> x, IList<double> a)
		{
			return SoftmaxOp_KM11_Sparse.A(x, a);
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="a">Buffer 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [Stochastic, SkipIfAllUniform] IList<Gaussian> x, IList<double> a, GaussianList result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			IList<Gaussian> prior = x.IsSparse() ? (IList<Gaussian>)SparseList<Gaussian>.FromSize(K) : (IList<Gaussian>)(new List<Gaussian>());
			prior.SetTo(new Func<Gaussian, Gaussian, Gaussian>((p, q) => p / q).Map(x, result));
			int evalCounter = 0;
			// z[1..K]: posterior means, z[K+1..2K]: log posterior variances
			var z = new Vector[2];
			z[0] = Vector.Zero(K, x.IsSparse() ? Sparsity.Sparse : Sparsity.Dense);
			z[0].SetTo(x.Select(o => o.GetMean()));
			z[1] = Vector.Zero(K, x.IsSparse() ? Sparsity.Sparse : Sparsity.Dense);
			z[1].SetTo(x.Select(o => Math.Log(o.GetVariance())));
			double startingValue = GradientAndValueAtPoint(prior, a, z[0], z[1], softmax.PseudoCount, null);
			var s = new LBFGSArray(5);
			s.MaximumStep = 1e3;
			s.MaximumIterations = 100;
			s.Epsilon = 1e-10;
			s.convergenceCriteria = BFGS.ConvergenceCriteria.Objective;
			z = s.Run(z, 1.0, delegate(Vector[] y, ref Vector[] grad) { evalCounter++; return GradientAndValueAtPoint(prior, a, y[0], y[1], softmax.PseudoCount, grad); });
			Func<double, double, Gaussian, Gaussian> f = (xmk, lxvk, priork) => {
				if (!priork.IsPointMass)
					return Gaussian.FromMeanAndVariance(xmk, Math.Exp(lxvk)) / priork;
				else
					return new Gaussian();
			};
			result.SetTo(f.Map(z[0], z[1], prior));
			double endValue = GradientAndValueAtPoint(prior, a, z[0], z[1], softmax.PseudoCount, null);
			//Console.WriteLine("Went from {0} to {1} in {2} steps, {3} evals", startingValue, endValue, s.IterationsPerformed, evalCounter);
			if (startingValue < endValue)
				Console.WriteLine("Warning: LBFGS resulted in an increased objective function");
			return result;
		}
	}


	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// Here the marginal prototype for logOdds is Gamma, which allows for heavier tailed distributions. 
	/// This implementation uses the Taylor series bound from Blei 06, followed by 
	/// optimization using LBFGS. This approach is linear in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Experimental)]
	public static class GammaSoftmaxOp
	{
		///<summary>
		///Helper function to get the means and variances of a list of Gaussians
		///</summary>
		private static void GetShapeAndRate(IList<Gamma> x, ref double[] a, ref double[] b)
		{
			int K = a.Length;
			for (int k = 0; k < K; k++) {
				a[k] = x[k].Shape;
				b[k] = x[k].Rate;
			}
		}


		//<summary>
		//Function to evaluate this factor's KL divergence contribution, and gradient (if grad is not null). 
		//</summary>
		//<param name="mu">Prior means</param>
		//<param name="s2">Prior variances</param>
		//<param name="x">x[1..K]: posterior mean, x[K+1..2K]: posterior log variance</param>
		//<param name="ns">Dirichlet counts-1</param>
		//<param name="grad">Vector to store the gradient in</param>
		//<returns>KL divergence</returns>
		private static double GradientAndValueAtPoint(double[] a2, double[] b2, Vector x, Vector ns, Vector grad)
		{
			int K = x.Count / 2;
			var a = new double[K];
			var b = new double[K];
			var exp = new double[K];
			for (int k = 0; k < K; k++) {
				a[k] = Math.Exp(x[k]);
				//might be better to use b=exp(m)+1
				b[k] = Math.Exp(x[k + K]) + 1;
				exp[k] = a[k] * (Math.Log(b[k]) - Math.Log(b[k] - 1));
			}
			double sum_n = ns.Sum();

			double logSumExp = MMath.LogSumExp(exp);
			double kl_value = sum_n * logSumExp;
			for (int k = 0; k < K; k++) {

				kl_value += Math.Log(b[k]) - MMath.GammaLn(a[k]) + (a[k] - 1) * MMath.Digamma(a[k]) - a[k] // entropy
                    - (a2[k] * Math.Log(b2[k]) - MMath.GammaLn(a2[k]) + (a2[k] - 1) * (MMath.Digamma(a[k]) - Math.Log(b[k])) - b2[k] * a[k] / b[k]) // cross entropy
                    - ns[k] * a[k] / b[k]; // factor
				if (double.IsInfinity(kl_value) || double.IsNaN(kl_value))
					throw new ApplicationException("doh");
				if (grad != null) {
					grad[k] = (a[k] - 1.0) * MMath.Trigamma(a[k]) - 1.0 // entropy
                        - (a2[k] - 1) * MMath.Trigamma(a[k]) + b2[k] / b[k] // cross term
                       - ns[k] / b[k] + sum_n * Math.Exp(exp[k] - logSumExp) * (Math.Log(b[k]) - Math.Log(b[k] - 1)); // factor
					grad[k] *= a[k]; // chain rule
					grad[K + k] = 1.0 / b[k]
                        + (a2[k] - 1) / b[k] - b2[k] * a[k] / (b[k] * b[k])
                       + ns[k] * a[k] / (b[k] * b[k]) + sum_n * Math.Exp(exp[k] - logSumExp) * a[k] * (1.0 / b[k] - 1.0 / (b[k] - 1));
					grad[K + k] *= b[k]; // chain rule
				}
			}

			return kl_value;
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static Vector GetMeanLog(IList<Gamma> x)
		{
			int K = x.Count;
			var a = new double[K];
			var b = new double[K];
			var m = new double[K];
			GetShapeAndRate(x, ref a, ref b);
			var exp = new double[K];
			for (int k = 0; k < K; k++) {
				exp[k] = a[k] * (Math.Log(b[k]) - Math.Log(b[k] - 1));
				m[k] = a[k] / b[k]; // means
			}
			double sum = MMath.LogSumExp(exp);
			var result = Vector.Constant(K, -sum);
			result += Vector.FromArray(m);
			return result;
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'softmax' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm([SkipIfAllUniform] IList<Gamma> x)
		{
			return Dirichlet.FromMeanLog(GetMeanLog(x));
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GammaList XAverageLogarithm<GammaList>([SkipIfUniform] Dirichlet softmax, [SkipIfAllUniform] GammaList x, GammaList result)
						 where GammaList : IList<Gamma>
		{
			int K = x.Count;
			var a2 = new double[K];
			var b2 = new double[K];
			var prior = new Gamma[K];
			for (int k = 0; k < K; k++)
				prior[k] = x[k] / result[k];
			GetShapeAndRate(prior, ref a2, ref b2);
			var a = new double[K];
			var b = new double[K];
			GetShapeAndRate(x, ref a, ref b);
			var counts = softmax.PseudoCount - 1;
			int evalCounter = 0;
			//z[1..K]: posterior means, z[K+1..2K]: log posterior variances
			var z = Vector.Zero(K * 2);
			for (int k = 0; k < K; k++) {
				z[k] = Math.Log(a[k]);
				z[K + k] = Math.Log(b[k] - 1);
			}
			double startingValue = GradientAndValueAtPoint(a2, b2, z, counts, null);
			var s = new LBFGS(5);
			s.debug = true;
			s.MaximumStep = 10;
			s.MaximumIterations = 100;
			s.Epsilon = 1e-5;
			s.convergenceCriteria = BFGS.ConvergenceCriteria.Objective;
			FunctionEval f = delegate(Vector y, ref Vector grad) { evalCounter++; return GradientAndValueAtPoint(a2, b2, y, counts, grad); };
			//DerivativeChecker.CheckDerivatives(f, z); 
			z = s.Run(z, 1.0, f);
			for (int k = 0; k < K; k++) {
				a[k] = Math.Exp(z[k]);
				b[k] = Math.Exp(z[K + k]) + 1;
				var rk = result[k];
				rk.SetToRatio(Gamma.FromShapeAndRate(a[k], b[k]), prior[k]);
				result[k] = rk;
			}
			double endValue = GradientAndValueAtPoint(a2, b2, z, counts, null);
			//Console.WriteLine("Went from {0} to {1} in {2} steps, {3} evals", startingValue, endValue, s.IterationsPerformed, evalCounter);
			if (startingValue < endValue)
				Console.WriteLine("Warning: LBFGS resulted in an increased objective function");
			return result;
		}
	}


	//[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	//public static class ExponentialSoftmaxOp
	//{
	//    /// <summary>
	//    /// Helper function to get the means and variances of a list of Gaussians
	//    /// </summary>
	//    private static void GetShapeAndRate<GammaList>(GammaList x, IList<double> a, IList<double> b)
	//        where GammaList : IList<Gamma>
	//    {
	//        int K = a.Count;
	//        for (int k = 0; k < K; k++)
	//        {
	//            a[k] = x[k].Shape;
	//            b[k] = x[k].Rate;
	//        }
	//    }

	//    /// <summary>
	//    /// Function to evaluate this factor's KL divergence contribution, and gradient (if grad is not null). 
	//    /// </summary>
	//    /// <param name="mu">Prior means</param>
	//    /// <param name="s2">Prior variances</param>
	//    /// <param name="x">x[1..K]: posterior mean, x[K+1..2K]: posterior log variance</param>
	//    /// <param name="ns">Dirichlet counts-1</param>
	//    /// <param name="grad">Vector to store the gradient in</param>
	//    /// <returns>KL divergence</returns>
	//    private static double GradientAndValueAtPoint(double[] b2, IList<double> x, Vector ns, Vector grad)
	//    {
	//        int K = x.Count;
	//        var b = new double[K];
	//        var exp = new double[K];
	//        for (int k = 0; k < K; k++)
	//        {
	//            // might be better to use b=exp(m)+1
	//            b[k] = Math.Exp(x[k]) + 1;
	//            exp[k] = Math.Log(b[k]) - Math.Log(b[k] - 1);
	//        }
	//        double sum_n = ns.Sum();

	//        double logSumExp = MMath.LogSumExp(exp);
	//        double kl_value = sum_n * logSumExp;
	//        for (int k = 0; k < K; k++)
	//        {

	//            kl_value += Math.Log(b[k]) - 1 // entropy
	//                - (Math.Log(b2[k]) - b2[k] / b[k]) // cross entropy
	//                - ns[k] / b[k]; // factor
	//            if (double.IsInfinity(kl_value) || double.IsNaN(kl_value))
	//                throw new ApplicationException("doh");
	//            if (grad != null)
	//            {
	//                grad[k] = 1.0 / b[k]
	//                   - b2[k] / (b[k] * b[k])
	//                   + ns[k] / (b[k] * b[k]) + sum_n * Math.Exp(exp[k] - logSumExp) * (1.0 / b[k] - 1.0 / (b[k] - 1));
	//                grad[k] *= (b[k] - 1.0); // chain rule
	//            }
	//        }

	//        return kl_value;
	//    }

	//    public static Vector GetMeanLog(IList<Gamma> x)
	//    {
	//        int K = x.Count;
	//        var a = new double[K];
	//        var b = new double[K];
	//        var m = new double[K];
	//        GetShapeAndRate(x, a, b);
	//        var exp = new double[K];
	//        for (int k = 0; k < K; k++)
	//        {
	//            exp[k] = Math.Log(b[k]) - Math.Log(b[k] - 1);
	//            m[k] = 1 / b[k]; // means
	//        }
	//        double sum = MMath.LogSumExp(exp);
	//        var result = Vector.Constant(K, -sum);
	//        result += Vector.FromArray(m);
	//        return result;
	//    }

	//    /// <summary>
	//    /// Message to softmax child node calculated using the bound used in Blei 2006
	//    /// E[log sum_i exp(x)] <= log E[sum_i exp(x)] = log sum_i exp(m[i]+v[i]/2)
	//    /// </summary>
	//    /// <param name="x">Incoming message from the logOdds, x</param>
	//    /// <returns></returns>
	//    public static Dirichlet SoftmaxAverageLogarithm([SkipIfAllUniform] IList<Gamma> x)
	//    {
	//        return Dirichlet.FromMeanLog(GetMeanLog(x));
	//    }

	//    /// <summary>
	//    /// VMP message to X from the softmax factor. The KL divergence is upper bounded using the bound used in Blei 2006:
	//    /// E[log sum_i exp(x)] <= log E[sum_i exp(x)] = log sum_i exp(m[i]+v[i]/2)
	//    /// LBFGS is used to jointly minimise over the posterior means and variances
	//    /// </summary>
	//    /// <param name="softmax">Incoming message from child</param>
	//    /// <param name="x">Current variational posterior for x</param>
	//    /// <param name="result">Previous output message</param>
	//    /// <returns>VMP message to X</returns>
	//    public static GammaList XAverageLogarithm<GammaList>([SkipIfUniform] Dirichlet softmax, [SkipIfAllUniform] GammaList x, GammaList result)
	//         where GammaList : IList<Gamma>
	//    {
	//        int K = x.Count;
	//        // a's should always be 1
	//        var a2 = new double[K];
	//        var b2 = new double[K];
	//        var prior = new Gamma[K];
	//        for (int k = 0; k < K; k++)
	//            prior[k] = x[k] / result[k];
	//        GetShapeAndRate(prior, a2, b2);
	//        var a = Vector.Zero(K);
	//        var b = Vector.Zero(K);
	//        GetShapeAndRate(x, a, b);
	//        var counts = softmax.PseudoCount - 1;
	//        int evalCounter = 0;
	//        // z[1..K]: posterior means, z[K+1..2K]: log posterior variances
	//        var z = Vector.Zero(K);
	//        z.SetToFunction(b, q => Math.Log(q - 1));
	//        double startingValue = GradientAndValueAtPoint(b2, z, counts, null);
	//        var s = new LBFGS(5);
	//        s.debug = true;
	//        s.MaximumStep = 10;
	//        s.MaximumIterations = 100;
	//        s.Epsilon = 1e-5;
	//        s.convergenceCriteria = BFGS.ConvergenceCriteria.Objective;
	//        FunctionEval f = delegate(Vector y, ref Vector grad) { evalCounter++; return GradientAndValueAtPoint(b2, y, counts, grad); };
	//        //DerivativeChecker.CheckDerivatives(f, z); 
	//        z = s.Run(z, 1.0, f);
	//        b.SetToFunction(z, q => Math.Exp(q) + 1);
	//        for (int k = 0; k < K; k++)
	//        {
	//            var rk = result[k];
	//            rk.SetToRatio(Gamma.FromShapeAndRate(1, b[k]), prior[k]);
	//            result[k] = rk;
	//        }
	//        double endValue = GradientAndValueAtPoint(b2, z, counts, null);
	//        Console.WriteLine("Went from {0} to {1} in {2} steps, {3} evals", startingValue, endValue, s.IterationsPerformed, evalCounter);
	//        if (startingValue < endValue)
	//            Console.WriteLine("Warning: LBFGS resulted in an increased objective function");
	//        return result;
	//    }
	//}

	/// <summary>
	/// Provides outgoing messages for <see cref="MMath.Softmax(IList{double})"/>, given random arguments to the function.
	/// This implementation uses the bound
	/// <c>exp(x[k])/sum_i exp(x[i]) >= prod_{i!=k} sigma(x[k]-x[i])</c>
	/// noted by Tom Minka. This approach is quadratic in the dimension K. 
	/// </summary>
	[FactorMethod(typeof(MMath), "Softmax", typeof(IList<double>))]
	[Quality(QualityBand.Preview)]
	public static class SoftmaxOp_ProductOfLogistics
	{
		//-- VMP -------------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(softmax,x))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor<GaussianList>([SkipIfAllUniform] GaussianList x, Dirichlet softmax, Dirichlet to_softmax)
				where GaussianList : IList<Gaussian>
		{
			double res = 0;
			for (int k = 0; k < softmax.Dimension; k++) {
				double z = 0;
				double mk, vk;
				x[k].GetMeanAndVariance(out mk, out vk);
				for (int j = 0; j < softmax.Dimension; j++) {
					if (j == k) continue;
					double mj, vj;
					x[j].GetMeanAndVariance(out mj, out vj);
					Gaussian logOdds = Gaussian.FromMeanAndVariance(mk - mj, vk + vj);
					z += BernoulliFromLogOddsOp_JJ96.AverageLogFactor(true, logOdds);
				}
				res += (softmax.PseudoCount[k] - 1.0) * z;
			}
			return res -softmax.GetLogNormalizer() - to_softmax.GetAverageLog(softmax);
		}

		/// <summary>
		/// VMP message to 'softmax'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'softmax' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(softmax,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Dirichlet SoftmaxAverageLogarithm<GaussianList>([SkipIfAllUniform] GaussianList x, Dirichlet result)
				where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			Vector meanLog = Vector.Zero(K);
			for (int k = 0; k < K; k++) {
				double z = 0;
				double mk, vk;
				x[k].GetMeanAndVariance(out mk, out vk);
				for (int j = 0; j < K; j++) {
					if (j == k) continue;
					double mj, vj;
					x[j].GetMeanAndVariance(out mj, out vj);
					Gaussian logOdds = Gaussian.FromMeanAndVariance(mk - mj, vk + vj);
					z += BernoulliFromLogOddsOp_JJ96.AverageLogFactor(true, logOdds);
				}
				meanLog[k] = z;
			}
			result.SetMeanLog(meanLog);
			return result;
		}

		const string NotSupportedMessage = "Variational Message Passing does not support a softmax factor with fixed output.";
		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'.</param>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static GaussianList XAverageLogarithm<GaussianList>(IList<double> softmax, IList<Gaussian> x, GaussianList result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="softmax">Incoming message from 'softmax'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'softmax' integrated out.
		/// The formula is <c>sum_softmax p(softmax) factor(softmax,x)</c>.
		/// </para>
		/// <para>The bound used below is exp(x[k])/sum_i exp(x[i]) >= prod_{i!=k} sigma(x[k]-x[i]).
		/// If UseBlei06Bound is set to true the bound E[log sum_i exp(x)] &lt;= log E[sum_i exp(x)]
		/// is used instead, followed by numerical optimisation using the LBFGS algorithm. </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="softmax"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static GaussianList XAverageLogarithm<GaussianList>([SkipIfUniform] Dirichlet softmax, [SkipIfAllUniform] GaussianList x, GaussianList result)
						where GaussianList : IList<Gaussian>
		{
			int K = x.Count;
			for (int k = 0; k < K; k++) {
				Gaussian gk = result[k];
				gk.SetToUniform();
				result[k] = gk;
			}
			for (int k = 0; k < K; k++) {
				Gaussian gk = result[k];
				double mk, vk;
				x[k].GetMeanAndVariance(out mk, out vk);
				double ak1 = softmax.PseudoCount[k] - 1;
				for (int j = 0; j < K; j++) {
					if (j == k) continue;
					Gaussian gj = result[j];
					double mj, vj;
					x[j].GetMeanAndVariance(out mj, out vj);
					Gaussian logProb = Gaussian.FromMeanAndVariance(mk - mj, vk + vj);
					Gaussian toLogProb;
					if (true) {
						toLogProb = BernoulliFromLogOddsOp_JJ96.LogOddsAverageLogarithm(true, logProb);
					} else {
						toLogProb = logProb;
						// iterate to convergence
						for (int iter = 0; iter < 20; iter++) {
							toLogProb = BernoulliFromLogOddsOp.LogOddsAverageLogarithm(true, logProb, toLogProb);
						}
					}
					if (false) {
						gk.Precision += ak1 * toLogProb.Precision;
						gk.MeanTimesPrecision += ak1 * toLogProb.MeanTimesPrecision;
						gj.Precision += ak1 * toLogProb.Precision;
						gj.MeanTimesPrecision -= ak1 * toLogProb.MeanTimesPrecision;
					} else {
						gk.SetToProduct(gk, DoubleMinusVmpOp.AAverageLogarithm(toLogProb, x[j])^ak1);
						if (double.IsNaN(gk.MeanTimesPrecision)) throw new Exception();
						gj.SetToProduct(gj, DoubleMinusVmpOp.BAverageLogarithm(toLogProb, x[k])^ak1);
					}
					result[j] = gj;
				}
				result[k] = gk;
			}
			return result;
		}
	}
}
