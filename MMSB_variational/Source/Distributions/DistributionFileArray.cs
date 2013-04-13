#define SpecializeInterfaces
using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using MicrosoftResearch.Infer.Collections;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Factors;
using System.IO;

namespace MicrosoftResearch.Infer.Distributions
{
	[Serializable]
	[Quality(QualityBand.Experimental)]
	public class DistributionFileArray<T, DomainType> : FileArray<T>, IDistribution<DomainType[]>, Sampleable<DomainType[]>
#if SpecializeInterfaces
, SettableTo<DistributionFileArray<T, DomainType>>,
		SettableToProduct<DistributionFileArray<T, DomainType>>,
		SettableToRatio<DistributionFileArray<T, DomainType>>,
		SettableToPower<DistributionFileArray<T, DomainType>>,
		SettableToWeightedSum<DistributionFileArray<T, DomainType>>,
		CanGetLogAverageOf<DistributionFileArray<T, DomainType>>,
		CanGetLogAverageOfPower<DistributionFileArray<T, DomainType>>,
		CanGetAverageLog<DistributionFileArray<T, DomainType>>
#endif
 where T : ICloneable, SettableTo<T>,
		SettableToProduct<T>,
		SettableToRatio<T>,
		SettableToPower<T>,
		SettableToWeightedSum<T>,
		CanGetLogAverageOf<T>,
		CanGetLogAverageOfPower<T>,
		CanGetAverageLog<T>,
		IDistribution<DomainType>,
		Sampleable<DomainType>
	{
		public DistributionFileArray(string prefix, int count, [SkipIfUniform] Func<int, T> init)
			: base(prefix, count, init)
		{
		}
		[Skip]
		public DistributionFileArray(string prefix, int count)
			: base(prefix, count)
		{
		}
		[Skip]
		public DistributionFileArray([IgnoreDeclaration] FileArray<DistributionFileArray<T, DomainType>> parent, int index, int count)
			: this(parent.GetItemFolder(index), count)
		{
			parent.containsFileArrays = true;
			this.doNotDelete = true;
			parent.StoreItem(index, this);
		}

		public object Clone()
		{
			if (containsFileArrays) throw new NotImplementedException();
			string folder = prefix.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar);
			return new DistributionFileArray<T,DomainType>(folder+"_clone", count, i => this[i]);
		}

		public void SetTo(DistributionFileArray<T, DomainType> value)
		{
			for (int i = 0; i < count; i++) {
				this[i] = value[i];
			}
		}

		public void SetToProduct(DistributionFileArray<T, DomainType> a, DistributionFileArray<T, DomainType> b)
		{
			for (int i = 0; i < count; i++) {
				T item = this[i];
				item.SetToProduct(a[i], b[i]);
				this[i] = item;
			}
		}

		public void SetToRatio(DistributionFileArray<T, DomainType> numerator, DistributionFileArray<T, DomainType> denominator)
		{
			throw new NotImplementedException();
		}

		public void SetToPower(DistributionFileArray<T, DomainType> value, double exponent)
		{
			throw new NotImplementedException();
		}

		public void SetToSum(double weight1, DistributionFileArray<T, DomainType> value1, double weight2, DistributionFileArray<T, DomainType> value2)
		{
			throw new NotImplementedException();
		}

		public double GetLogAverageOf(DistributionFileArray<T, DomainType> that)
		{
			throw new NotImplementedException();
		}

		public double GetLogAverageOfPower(DistributionFileArray<T, DomainType> that, double power)
		{
			throw new NotImplementedException();
		}

		public double GetAverageLog(DistributionFileArray<T, DomainType> that)
		{
			throw new NotImplementedException();
		}

		public void SetToUniform()
		{
			for (int i = 0; i < count; i++) {
				T item = this[i];
				item.SetToUniform();
				this[i] = item;
			}
		}

		public bool IsUniform()
		{
			return this.All(item => item.IsUniform());
		}

		public DomainType[] Point
		{
			get
			{
				throw new NotImplementedException();
			}
			set
			{
				throw new NotImplementedException();
			}
		}

		public bool IsPointMass
		{
			get { throw new NotImplementedException(); }
		}

		public double MaxDiff(object that)
		{
			throw new NotImplementedException();
		}

		public double GetLogProb(DomainType[] value)
		{
			throw new NotImplementedException();
		}

		public DomainType[] Sample()
		{
			throw new NotImplementedException();
		}

		public DomainType[] Sample(DomainType[] result)
		{
			throw new NotImplementedException();
		}
	}
}
