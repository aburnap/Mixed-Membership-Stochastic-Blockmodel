// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides factors and operators for using Enum types.
	/// </summary>
	public class EnumSupport
	{

		/// <summary>
		/// Converts an Enum to an Int
		/// </summary>
		/// <param name="en"></param>
		/// <returns></returns>
		[ParameterNames("Int", "Enum")]
		public static int EnumToInt<TEnum>(TEnum en)
		{
			return (int)(object)en;
		}

		/// <summary>
		/// Samples an enum value from a discrete distribution.
		/// </summary>
		/// <typeparam name="TEnum">The type of the enum to sample</typeparam>
		/// <param name="probs">Vector of the probability of each Enum value, in order</param>
		/// <returns>An enum sampled from the distribution</returns>
		[Stochastic]
		[ParameterNames("Sample", "Probs")]
		public static TEnum DiscreteEnum<TEnum>(Vector probs)
		{
			int i= MicrosoftResearch.Infer.Distributions.Discrete.Sample(probs);
			return (TEnum)Enum.GetValues(typeof(TEnum)).GetValue(i);
		}

		/// <summary>
		/// Test if two enums are equal.
		/// </summary>
		/// <param name="a">First integer</param>
		/// <param name="b">Second integer</param>
		/// <returns>True if a==b.</returns>
		public static bool AreEqual<TEnum>(TEnum a, TEnum b) { return EnumToInt<TEnum>(a) == EnumToInt<TEnum>(b); }
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="EnumSupport.EnumToInt{TEnum}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(EnumSupport), "EnumToInt<>")]
	[Quality(QualityBand.Preview)]
	public static class EnumToIntOp<TEnum>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="Enum">Constant value for 'Enum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Int,Enum))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int Int, TEnum Enum)
		{
			return (EnumSupport.EnumToInt(Enum) == Int) ? 0.0 : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="Enum">Constant value for 'Enum'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Int,Enum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int Int, TEnum Enum) { return LogAverageFactor(Int, Enum); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="Enum">Constant value for 'Enum'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Int,Enum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(int Int, TEnum Enum) { return LogAverageFactor(Int, Enum); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="Enum">Incoming message from 'Enum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Enum) p(Enum) factor(Int,Enum))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int Int, DiscreteEnum<TEnum> Enum)
		{
			return Enum.GetLogProb(Int);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="Enum">Incoming message from 'Enum'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Enum) p(Enum) factor(Int,Enum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int Int, DiscreteEnum<TEnum> Enum)
		{
			return LogAverageFactor(Int, Enum);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Incoming message from 'Int'.</param>
		/// <param name="Enum">Constant value for 'Enum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Int) p(Int) factor(Int,Enum))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete Int, TEnum Enum)
		{
			return Int.GetLogProb(EnumSupport.EnumToInt(Enum));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Incoming message from 'Int'.</param>
		/// <param name="Enum">Incoming message from 'Enum'.</param>
		/// <param name="to_Int">Outgoing message to 'Int'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Int,Enum) p(Int,Enum) factor(Int,Enum))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete Int, DiscreteEnum<TEnum> Enum, [Fresh] Discrete to_Int)
		{
			return to_Int.GetLogAverageOf(Int);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Int">Incoming message from 'Int'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Int) p(Int) factor(Int,Enum) / sum_Int p(Int) messageTo(Int))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Discrete Int) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Int,Enum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// EP message to 'Int'
		/// </summary>
		/// <param name="Enum">Incoming message from 'Enum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Int' as the random arguments are varied.
		/// The formula is <c>proj[p(Int) sum_(Enum) p(Enum) factor(Int,Enum)]/p(Int)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Enum"/> is not a proper distribution</exception>
		public static Discrete IntAverageConditional([SkipIfUniform] DiscreteEnum<TEnum> Enum, Discrete result)
		{
			result.SetTo(Enum);
			return result;
		}
		[Skip]
		public static Discrete IntAverageConditionalInit([IgnoreDependency] DiscreteEnum<TEnum> Enum)
		{
			return Discrete.Uniform(Enum.Dimension);
		}

		/// <summary>
		/// EP message to 'Enum'
		/// </summary>
		/// <param name="Int">Incoming message from 'Int'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Enum' as the random arguments are varied.
		/// The formula is <c>proj[p(Enum) sum_(Int) p(Int) factor(Int,Enum)]/p(Enum)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Int"/> is not a proper distribution</exception>
		public static DiscreteEnum<TEnum> EnumAverageConditional([SkipIfUniform] Discrete Int, DiscreteEnum<TEnum> result)
		{
			result.SetTo(Int);
			return result;
		}

		/// <summary>
		/// EP message to 'Enum'
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Enum' conditioned on the given values.
		/// </para></remarks>
		public static DiscreteEnum<TEnum> EnumAverageConditional(int Int, DiscreteEnum<TEnum> result)
		{
			result.Point = Int;
			return result;
		}

		/// <summary>
		/// VMP message to 'Int'
		/// </summary>
		/// <param name="Enum">Incoming message from 'Enum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Int' as the random arguments are varied.
		/// The formula is <c>proj[sum_(Enum) p(Enum) factor(Int,Enum)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Enum"/> is not a proper distribution</exception>
		public static Discrete IntAverageLogarithm([SkipIfUniform] DiscreteEnum<TEnum> Enum, Discrete result)
		{
			return IntAverageConditional(Enum, result);
		}

		/// <summary>
		/// VMP message to 'Enum'
		/// </summary>
		/// <param name="Int">Incoming message from 'Int'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Enum' with 'Int' integrated out.
		/// The formula is <c>sum_Int p(Int) factor(Int,Enum)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Int"/> is not a proper distribution</exception>
		public static DiscreteEnum<TEnum> EnumAverageLogarithm([SkipIfUniform] Discrete Int, DiscreteEnum<TEnum> result)
		{
			return EnumAverageConditional(Int, result);
		}

		/// <summary>
		/// VMP message to 'Enum'
		/// </summary>
		/// <param name="Int">Constant value for 'Int'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Enum' conditioned on the given values.
		/// </para></remarks>
		public static DiscreteEnum<TEnum> EnumAverageLogarithm(int Int, DiscreteEnum<TEnum> result)
		{
			return EnumAverageConditional(Int, result);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="EnumSupport.DiscreteEnum&lt;TEnum&gt;"/>, given random arguments to the function.
	/// </summary>
	/// <remarks>
	/// This class provides operators which have Enum arguments.  
	/// The rest are provided by DiscreteFromDirichlet.
	/// </remarks>
	[FactorMethod(typeof(EnumSupport), "DiscreteEnum<>")]
	[Quality(QualityBand.Stable)]
	public static class DiscreteEnumFromDirichletOp<TEnum>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="probs">Incoming message from 'Probs'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Probs) p(Probs) factor(Sample,Probs))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(TEnum sample, Dirichlet probs)
		{
			return DiscreteFromDirichletOp.LogAverageFactor(EnumSupport.EnumToInt(sample), probs);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="probs">Constant value for 'Probs'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sample,Probs))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(TEnum sample, Vector probs)
		{
			return DiscreteFromDirichletOp.LogAverageFactor(EnumSupport.EnumToInt(sample), probs);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="probs">Incoming message from 'Probs'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(Probs) p(Probs) log(factor(Sample,Probs))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(TEnum sample, Dirichlet probs)
		{
			return DiscreteFromDirichletOp.AverageLogFactor(EnumSupport.EnumToInt(sample), probs);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="probs">Constant value for 'Probs'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sample,Probs))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(TEnum sample, Vector probs)
		{
			return DiscreteFromDirichletOp.AverageLogFactor(EnumSupport.EnumToInt(sample), probs);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="probs">Incoming message from 'Probs'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Probs) p(Probs) factor(Sample,Probs))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(TEnum sample, Dirichlet probs)
		{
			return DiscreteFromDirichletOp.LogEvidenceRatio(EnumSupport.EnumToInt(sample), probs);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="probs">Constant value for 'Probs'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sample,Probs))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(TEnum sample, Vector probs)
		{
			return DiscreteFromDirichletOp.LogEvidenceRatio(EnumSupport.EnumToInt(sample), probs);
		}

		/// <summary>
		/// EP message to 'Probs'
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Probs' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet ProbsAverageConditional(TEnum sample, Dirichlet result)
		{
			return DiscreteFromDirichletOp.ProbsAverageConditional(EnumSupport.EnumToInt(sample), result);
		}
		/// <summary>
		/// VMP message to 'Probs'
		/// </summary>
		/// <param name="sample">Constant value for 'Sample'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Probs' conditioned on the given values.
		/// </para></remarks>
		public static Dirichlet ProbsAverageLogarithm(TEnum sample, Dirichlet result)
		{
			return DiscreteFromDirichletOp.ProbsAverageLogarithm(EnumSupport.EnumToInt(sample), result);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="EnumSupport.AreEqual{TEnum}"/>, given random arguments to the function.
	/// </summary>
	/// <remarks>This class only implements enum-specific overloads that are not provided by DiscreteAreEqualOp<para></para></remarks>
	[FactorMethod(typeof(EnumSupport), "AreEqual<>")]
	[Quality(QualityBand.Stable)]
	public class DiscreteEnumAreEqualOp<TEnum>
	{
		private static int ToInt(TEnum en) { return (int)(object)en; }

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(areEqual) p(areEqual) factor(areEqual,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static DiscreteEnum<TEnum> AAverageConditional([SkipIfUniform] Bernoulli areEqual, TEnum B, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.AAverageConditional(areEqual, ToInt(B), result));
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static DiscreteEnum<TEnum> AAverageConditional(bool areEqual, TEnum B, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.AAverageConditional(areEqual, ToInt(B), result));
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static DiscreteEnum<TEnum> AAverageLogarithm(bool areEqual, TEnum B, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.AAverageLogarithm(areEqual, ToInt(B), result));
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'areEqual' integrated out.
		/// The formula is <c>sum_areEqual p(areEqual) factor(areEqual,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static DiscreteEnum<TEnum> AAverageLogarithm([SkipIfUniform] Bernoulli areEqual, TEnum B, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.AAverageLogarithm(areEqual, ToInt(B), result));
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(areEqual) p(areEqual) factor(areEqual,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static DiscreteEnum<TEnum> BAverageConditional([SkipIfUniform] Bernoulli areEqual, TEnum A, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.BAverageConditional(areEqual, ToInt(A), result));
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static DiscreteEnum<TEnum> BAverageConditional(bool areEqual, TEnum A, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.BAverageConditional(areEqual, ToInt(A), result));
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'areEqual' integrated out.
		/// The formula is <c>sum_areEqual p(areEqual) factor(areEqual,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static DiscreteEnum<TEnum> BAverageLogarithm([SkipIfUniform] Bernoulli areEqual, TEnum A, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.BAverageLogarithm(areEqual, ToInt(A), result));
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static DiscreteEnum<TEnum> BAverageLogarithm(bool areEqual, TEnum A, Discrete result)
		{
			return new DiscreteEnum<TEnum>(DiscreteAreEqualOp.BAverageLogarithm(areEqual, ToInt(A), result));
		}
		/// <summary>
		/// EP message to 'areEqual'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[p(areEqual) sum_(b) p(b) factor(areEqual,a,b)]/p(areEqual)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageConditional(TEnum A, Discrete B)
		{
			return DiscreteAreEqualOp.AreEqualAverageConditional(ToInt(A), B);
		}
		/// <summary>
		/// EP message to 'areEqual'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[p(areEqual) sum_(a) p(a) factor(areEqual,a,b)]/p(areEqual)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageConditional(Discrete A, TEnum B)
		{
			return DiscreteAreEqualOp.AreEqualAverageConditional(A, ToInt(B));
		}
		/// <summary>
		/// VMP message to 'areEqual'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(areEqual,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageLogarithm(TEnum A, Discrete B)
		{
			return DiscreteAreEqualOp.AreEqualAverageLogarithm(ToInt(A), B);
		}
		/// <summary>
		/// VMP message to 'areEqual'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(areEqual,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageLogarithm(Discrete A, TEnum B)
		{
			return DiscreteAreEqualOp.AreEqualAverageLogarithm(A, ToInt(B));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual) p(areEqual) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli areEqual, TEnum A, TEnum B)
		{
			return DiscreteAreEqualOp.LogAverageFactor(areEqual, ToInt(A), ToInt(B));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, TEnum A, TEnum B)
		{
			return DiscreteAreEqualOp.LogAverageFactor(areEqual, ToInt(A), ToInt(B));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, TEnum A, Discrete B)
		{
			return DiscreteAreEqualOp.LogAverageFactor(areEqual, ToInt(A), B);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Discrete A, TEnum B)
		{
			return DiscreteAreEqualOp.LogAverageFactor(areEqual, A, ToInt(B));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_B">Outgoing message to 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, TEnum A, Discrete B, [Fresh] Discrete to_B)
		{
			return DiscreteAreEqualOp.LogEvidenceRatio(areEqual, ToInt(A), B, to_B);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="to_A">Outgoing message to 'A'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, Discrete A, TEnum B, [Fresh] Discrete to_A)
		{
			return DiscreteAreEqualOp.LogEvidenceRatio(areEqual, A, ToInt(B), to_A);
		}
		public static double LogEvidenceRatio(bool areEqual, TEnum A, TEnum B)
		{
			return DiscreteAreEqualOp.LogEvidenceRatio(areEqual, ToInt(A), ToInt(B));
		}
	}
}
