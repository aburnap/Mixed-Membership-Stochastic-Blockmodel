// (C) Copyright 2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Copy{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "Output", "Input" }, typeof(Factor), "OffsetCopy<>")]
	[Quality(QualityBand.Experimental)]
	public static class OffsetCopyOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(OffsetCopy,value))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="output">Incoming message from 'OffsetCopy'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(OffsetCopy) p(OffsetCopy) factor(OffsetCopy,value) / sum_OffsetCopy p(OffsetCopy) messageTo(OffsetCopy))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio<TDist>(TDist output)
			where TDist : IDistribution<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="output">Constant value for 'OffsetCopy'.</param>
		/// <param name="input">Incoming message from 'value'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(value) p(value) factor(OffsetCopy,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<TDist>(T output, TDist input)
			where TDist : CanGetLogProb<T>
		{
			return input.GetLogProb(output);
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="input">Incoming message from 'value'.</param>
		/// <param name="Output">Incoming message from 'OffsetCopy'.</param>
		/// <returns>The outgoing EP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(OffsetCopy) p(OffsetCopy) factor(OffsetCopy,value)]/p(value)</c>.
		/// </para></remarks>
		public static TDist InputAverageConditional<TDist, TInput>([RequiredArgument] TInput input, [IgnoreDependency] TDist Output)
		{
			return Output;
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="input">Incoming message from 'value'.</param>
		/// <param name="output">Constant value for 'OffsetCopy'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' conditioned on the given values.
		/// </para></remarks>
		public static TDist InputAverageConditional<TDist>([RequiredArgument] TDist input, [IgnoreDependency] T output, TDist result)
			where TDist : IDistribution<T>
		{
			result.Point = output;
			return result;
		}

		/// <summary>
		/// EP message to 'OffsetCopy'
		/// </summary>
		/// <param name="Output">Incoming message from 'OffsetCopy'.</param>
		/// <param name="Input">Incoming message from 'value'.</param>
		/// <returns>The outgoing EP message to the 'OffsetCopy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'OffsetCopy' as the random arguments are varied.
		/// The formula is <c>proj[p(OffsetCopy) sum_(value) p(value) factor(OffsetCopy,value)]/p(OffsetCopy)</c>.
		/// </para></remarks>
		public static TDist OutputAverageConditional<TDist, TOutput>([RequiredArgument] TOutput Output, [IgnoreDependency] TDist Input)
		{
			return Input;
		}

		//-- VMP ------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(OffsetCopy,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		// must have upward Trigger to match the Trigger on UsesEqualDef.UsesAverageLogarithm
		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="Input">Incoming message from 'value'.</param>
		/// <param name="Output">Incoming message from 'OffsetCopy'.</param>
		/// <returns>The outgoing VMP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' with 'OffsetCopy' integrated out.
		/// The formula is <c>sum_OffsetCopy p(OffsetCopy) factor(OffsetCopy,value)</c>.
		/// </para></remarks>
		public static TDist InputAverageLogarithm<TDist, TInput>([RequiredArgument] TInput Input, [IgnoreDependency]TDist Output)
		{
			return Output;
		}

		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="input">Incoming message from 'value'.</param>
		/// <param name="output">Constant value for 'OffsetCopy'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' conditioned on the given values.
		/// </para></remarks>
		public static TDist InputAverageLogarithm<TDist>([RequiredArgument] TDist input, [IgnoreDependency] T output, TDist result)
			where TDist : IDistribution<T>
		{
			result.Point = output;
			return result;
		}


		/// <summary>
		/// VMP message to 'OffsetCopy'
		/// </summary>
		/// <param name="Output">Incoming message from 'OffsetCopy'.</param>
		/// <param name="Input">Incoming message from 'value'.</param>
		/// <returns>The outgoing VMP message to the 'OffsetCopy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'OffsetCopy' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(OffsetCopy,value)]</c>.
		/// </para></remarks>
		public static TDist OutputAverageLogarithm<TDist, TOutput>([RequiredArgument] TOutput Output, [IgnoreDependency] TDist Input)
		{
			return Input;
		}

	}
}
