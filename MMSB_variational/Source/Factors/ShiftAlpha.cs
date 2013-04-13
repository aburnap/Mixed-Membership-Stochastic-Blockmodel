// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Factors that change a message channel's alpha factor
	/// </summary>
	[Hidden]
	public static class ShiftAlpha
	{
		/// <summary>
		/// Changes a message channel's alpha value, going to a factor.
		/// </summary>
		/// <typeparam name="T">The domain type.</typeparam>
		[ParameterNames("factor", "variable", "variableAlpha", "factorAlpha")]
		public static T ToFactor<T>(T variable, double variableAlpha, double factorAlpha) { return variable; }
		/// <summary>
		/// Changes a message channel's alpha value, coming from a factor.
		/// </summary>
		/// <typeparam name="T">The domain type.</typeparam>
		[ParameterNames("variable", "factor", "factorAlpha", "variableAlpha")]
		public static T FromFactor<T>(T factor, double factorAlpha, double variableAlpha) { return factor; }
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="ShiftAlpha.ToFactor{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(ShiftAlpha), "ToFactor<>")]
	[Quality(QualityBand.Experimental)]
	public static class ShiftAlphaToFactorOp
	{
		/// <summary>
		/// EP message to 'factor'
		/// </summary>
		/// <param name="factor">Incoming message from 'factor'.</param>
		/// <param name="variable">Incoming message from 'variable'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="variableAlpha">Constant value for 'variableAlpha'.</param>
		/// <param name="factorAlpha">Constant value for 'factorAlpha'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'factor' as the random arguments are varied.
		/// The formula is <c>proj[p(factor) sum_(variable) p(variable) factor(factor,variable,variableAlpha,factorAlpha)]/p(factor)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="variable"/> is not a proper distribution</exception>
		public static T FactorAverageConditional<T>(T factor, [SkipIfUniform] T variable, double variableAlpha, double factorAlpha, T result)
			where T : SettableToPower<T>, SettableToProduct<T>
		{
			result.SetToPower(factor, variableAlpha-factorAlpha);
			result.SetToProduct(result, variable);
			return result;
		}

		/// <summary>
		/// EP message to 'variable'
		/// </summary>
		/// <param name="factor">Incoming message from 'factor'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'variable' as the random arguments are varied.
		/// The formula is <c>proj[p(variable) sum_(factor) p(factor) factor(factor,variable,variableAlpha,factorAlpha)]/p(variable)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="factor"/> is not a proper distribution</exception>
		public static T VariableAverageConditional<T>([SkipIfUniform] T factor, T result)
			where T : SettableTo<T>
		{
			result.SetTo(factor);
			return result;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="factor">Incoming message from 'factor'.</param>
		/// <param name="variable">Incoming message from 'variable'.</param>
		/// <param name="variableAlpha">Constant value for 'variableAlpha'.</param>
		/// <param name="factorAlpha">Constant value for 'factorAlpha'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx / int ftilde(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx / int ftilde(x) qnotf(x) dx)</c>
		/// where <c>x = (factor,variable,variableAlpha,factorAlpha)</c>.
		/// </para></remarks>
		public static double LogEvidenceRatioOld<T>(T factor, T variable, double variableAlpha, double factorAlpha)
			where T : ICloneable, CanGetAverageLog<T>, SettableToPower<T>, SettableToProduct<T>
		{
			if (variableAlpha == 1 && factorAlpha == 0) {
				// EP variable to VMP factor
				T to_factor = FactorAverageConditional(factor, variable, variableAlpha, factorAlpha, (T)variable.Clone());
				return -to_factor.GetAverageLog(factor);
			} else if (variableAlpha == 0 && factorAlpha == 1) {
				// VMP variable to EP factor
				return variable.GetAverageLog(factor);
			} else {
				throw new NotImplementedException();
			}
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="ShiftAlpha.FromFactor{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(ShiftAlpha), "FromFactor<>")]
	[Quality(QualityBand.Experimental)]
	public static class ShiftAlphaFromFactorOp
	{
		/// <summary>
		/// EP message to 'factor'
		/// </summary>
		/// <param name="factor">Incoming message from 'factor'.</param>
		/// <param name="variable">Incoming message from 'variable'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="factorAlpha">Constant value for 'factorAlpha'.</param>
		/// <param name="variableAlpha">Constant value for 'variableAlpha'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'factor' as the random arguments are varied.
		/// The formula is <c>proj[p(factor) sum_(variable) p(variable) factor(variable,factor,factorAlpha,variableAlpha)]/p(factor)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="variable"/> is not a proper distribution</exception>
		public static T FactorAverageConditional<T>(T factor, [SkipIfUniform] T variable, double factorAlpha, double variableAlpha, T result)
			where T : SettableToPower<T>, SettableToProduct<T>
		{
			result.SetToPower(factor, variableAlpha-factorAlpha);
			result.SetToProduct(result, variable);
			return result;
		}

		/// <summary>
		/// EP message to 'variable'
		/// </summary>
		/// <param name="factor">Incoming message from 'factor'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'variable' as the random arguments are varied.
		/// The formula is <c>proj[p(variable) sum_(factor) p(factor) factor(variable,factor,factorAlpha,variableAlpha)]/p(variable)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="factor"/> is not a proper distribution</exception>
		public static T VariableAverageConditional<T>([SkipIfUniform] T factor, T result)
			where T : SettableTo<T>
		{
			result.SetTo(factor);
			return result;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="factor">Incoming message from 'factor'.</param>
		/// <param name="variable">Incoming message from 'variable'.</param>
		/// <param name="variableAlpha">Constant value for 'variableAlpha'.</param>
		/// <param name="factorAlpha">Constant value for 'factorAlpha'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx / int ftilde(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx / int ftilde(x) qnotf(x) dx)</c>
		/// where <c>x = (variable,factor,factorAlpha,variableAlpha)</c>.
		/// </para></remarks>
		public static double LogEvidenceRatioOld<T>(T factor, T variable, double variableAlpha, double factorAlpha)
			where T : ICloneable, CanGetAverageLog<T>, SettableToPower<T>, SettableToProduct<T>
		{
			if (variableAlpha == 1 && factorAlpha == 0) {
				// EP variable to VMP factor
				T to_factor = FactorAverageConditional(factor, variable, variableAlpha, factorAlpha, (T)variable.Clone());
				return -to_factor.GetAverageLog(factor);
			} else if (variableAlpha == 0 && factorAlpha == 1) {
				// VMP variable to EP factor
				return variable.GetAverageLog(factor);
			} else {
				throw new NotImplementedException();
			}
		}
	}
}
