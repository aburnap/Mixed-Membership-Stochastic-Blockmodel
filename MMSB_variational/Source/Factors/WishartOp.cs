// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Wishart.SampleFromShapeAndScale"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Wishart), "SampleFromShapeAndScale")]
	[Quality(QualityBand.Stable)]
	public static class WishartFromShapeAndScaleOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="scale">Constant value for 'scale'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,shape,scale))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(PositiveDefiniteMatrix sample, double shape, PositiveDefiniteMatrix scale)
		{
			Wishart to_sample = SampleAverageConditional(shape, scale);
			return to_sample.GetLogProb(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="scale">Constant value for 'scale'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,shape,scale))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(PositiveDefiniteMatrix sample, double shape, PositiveDefiniteMatrix scale) { return LogAverageFactor(sample, shape, scale); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="scale">Constant value for 'scale'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,shape,scale))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(PositiveDefiniteMatrix sample, double shape, PositiveDefiniteMatrix scale) { return LogAverageFactor(sample, shape, scale); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,shape,scale))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Wishart sample, [Fresh] Wishart to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="scale">Constant value for 'scale'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Wishart SampleAverageLogarithm(double shape, PositiveDefiniteMatrix scale)
		{
			return Wishart.FromShapeAndScale(shape, scale);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="scale">Constant value for 'scale'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Wishart SampleAverageConditional(double shape, PositiveDefiniteMatrix scale)
		{
			return Wishart.FromShapeAndScale(shape, scale);
		}
	}
}
