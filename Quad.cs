using System;
using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

/// <summary>
/// Copyright (c) Cristian Cozma 2024. Licensed under MIT License.
/// 
/// A light weight library for double double arithmetic.
/// It can incorporate any built in numeric data type, including decimal.
/// Unlike any other libraries, it uses a new class of EFT's :
/// Fused Error Free Transformations (FEFT).
/// </summary>

namespace Sabs.Numerics
{
    [Serializable, StructLayout(LayoutKind.Sequential), DebuggerDisplay("{ToStringInvariant()}")]
    public struct Quad : IEquatable<Quad>, IComparable<Quad>, IEquatable<double>, IComparable<double>
    {
        private double h, l;

        public static readonly Quad Zero = 0;
        public static readonly Quad NaN = double.NaN;
        public static readonly Quad Epsilon = double.Epsilon;
        public static readonly Quad PositiveInfinity = double.PositiveInfinity;
        public static readonly Quad NegativeInfinity = double.NegativeInfinity;
        public static readonly Quad MaxValue = new Quad() { h = double.MaxValue, l = double.MaxValue / 18014398509481984 };
        public static readonly Quad MinValue = new Quad() { h = double.MinValue, l = double.MinValue / 18014398509481984 };

        public static readonly Quad E; // = Parse("2.718281828459045235360287471352662498");
        // e1: 2.718281828459045, e2: 2.3536728122053305e-16, e3: -2.1643445851410196e-32, e4: -1.0478056262447993e-48
        public static readonly Quad PI; // = Parse("3.141592653589793238462643383279502884");
        // helper constant for argument range reduction for sin and cos
        private static readonly Quad pilow = new Quad() { h = -2.9947698097183397e-33, l = 1.1124542208633653e-49 };
        // pi1: 3.141592653589793, pi2: 1.2246467991473532e-16, pi3: -2.9947698097183397e-33, pi4: 1.1124542208633653e-49
        private static readonly Quad sqrt2; // = Parse(1.414213562373095048801688724209698");
        private static readonly Quad ln2; // = Parse("0.693147180559945309417232121458176568");
        // ln2_1: 0.6931471805599453, ln2_2: 2.3190468138462996e-17 ln2_3: 1.246116740312015e-33, ln2_4: 4.831093126581452e-50
        // helper constant 2/ln(2)
        private static readonly Quad inv2ln2; // = Parse("2.885390081777926814719849362003784");
        // helper constant 1/log2(10)
        private static readonly Quad invlg10; // = Parse("0.30102999566398119521373889472449");

        static Quad()
        {
            //E = new Quad(PowHelp(ln2 = 1), 1);
            E = Pow(PowHelp(ln2 = 0.5), 4);
            inv2ln2 = ~(ln2 = LogHelp2(sqrt2 = Sqrt(inv2ln2 = 2)));
            //inv2ln2.l -= 6.1629758220391547E-32;
            ln2.Scale(2);
            //PI = AtanHelp(~((Sqrt(3) + sqrt2) * (sqrt2 + 1))) * 24;
            PI = AtanHelp(~(Sqrt(4 + sqrt2 * 2) + (sqrt2 + 1))) * 16;
            invlg10 = ~Log2((Quad)10);
        }

        public Quad(string s) : this(s, NumberFormatInfo.InvariantInfo)
        { }

        public Quad(string s, NumberFormatInfo info)
        {
            if (s == null)
            {
                throw new ArgumentNullException("s");
            }
            if (!TryParse(s, info, out this))
            {
                if (IsFinite(h))
                {
                    throw new FormatException();
                }
                throw new ArgumentOutOfRangeException("s");
            }
        }

        // constructors for normalization
        public Quad(double vh, double vl)
        {
            h = vh + vl;
            l = (vh - h) + vl;
        }

        public Quad(Quad q, double v)
        {
            h = (v + q.l) + q.h;
            l = (q.h - h) + v + q.l;
        }

        public Quad(double v, Quad q)
        {
            h = (q.h + q.l) + v;
            l = (v - h) + q.h + q.l;
        }

        public static explicit operator double(Quad q)
        {
            return q.h + q.l;
        }

        public static explicit operator int(Quad q)
        {
            // need to go through long because int.MaxValue + 1 and int.MinValue - 1 fit in a double,
            // but both values may be valid as ints depending on the sign of the remainder
            return (int)q.ToLong();
        }

        private long ToLong()
        {
            long v = (long)(h + l);
            double t = (h - v) + l;
            if (v != 0 && t != 0)
            {
                return v - ((t < 0 ? 1 : 0) + (v >> 63));
            }
            return v;
        }

        public static explicit operator long(Quad q)
        {
            long v = (long)q.l;
            double t = q.l - v;
            if (q.h == long.MaxValue) // this is in fact long.MaxValue + 1
            {
                return v - ((t < 0 ? 1 : 0) + long.MinValue);
            }
            q.l = t;
            return v + q.ToLong();
        }

        public static decimal ToDecimal(Quad q, int decimals, MidpointRounding mode)
        {
            double v = Log2(q.h);
            long h, l;
            int n = 28;
            if (v > 0)
            {
                n -= (int)((v + 0.42198608378962976) * 0.3010299956639812);
            }
            if (n > decimals) n = decimals;

            int s = q.h < 0 ? -5 : 5;
            q.Scale(Pow2(n - 43));
            q *= Pow((double)s, n);

            q.h -= (h = (long)q.h);
            v = (q.h + q.l) * (1L << 43);

            // check for underestimate
            if (n < decimals && (h < 0x3333333333333L ||
                (h == 0x3333333333333L && v < 1759218604441.55)))
            {
                h *= 10;
                v *= 10;
                n++;
            }

            // round the last 4 bits
            //double t = (h - h * 0.9999999999999999) * 5;
            //v = (v + t) - t;
            l = (long)(mode == MidpointRounding.ToEven ? Math.Round(v) : (v + 0.5));
            h = (h << 11) + (l >> 32);

            return new decimal(unchecked((int)l), unchecked((int)h), (int)(h >> 32), s < 0, checked((byte)n));
        }

        public static explicit operator decimal(Quad q)
        {
            return ToDecimal(q, 28, MidpointRounding.ToEven);
        }

        public static implicit operator Quad(double v)
        {
            return new Quad { h = v };
        }

        public static implicit operator Quad(int v)
        {
            return new Quad { h = v };
        }

        public static implicit operator Quad(long v)
        {
            Quad r;
            r.h = v;
            if (r.h == long.MaxValue)
            {
                v += long.MinValue;
            }
            else
            {
                v -= (long)r.h;
            }
            r.l = v;
            return r;
        }

        public static implicit operator Quad(decimal v)
        {
            int[] i = decimal.GetBits(v);
            int n = i[3] >> 16;
            //uint[] b = (uint[])(object)i;

            return new Quad(unchecked(((long)(uint)i[2] << 21) + ((uint)i[1] >> 11)) * (double)(1L << 43),
            unchecked(((long)(i[1] & 0x7ff) << 32) + (uint)i[0])) * Pow(n < 0 ? -10 : 10, -(n & sbyte.MaxValue));
        }

        public override bool Equals(object obj)
        {
            return obj is Quad && Equals((Quad)obj);
        }

        public override int GetHashCode()
        {
            double v = h + l;
            return v.GetHashCode() ^ (h - v + l).GetHashCode();
        }

        public bool Equals(Quad x)
        {
            return h + l == x.h + x.l && l == (x.h - h) + x.l;
        }

        public bool Equals(double v)
        {
            return h + l == v && v - h == l;
        }

        public int CompareTo(Quad x)
        {
            double t = h + l;
            double v = x.h + x.l;
            if (t > v) return 1;
            if (t < v) return -1;
            v = (x.h - h) + x.l; t = l;
            if (t > v) return 1;
            if (t < v) return -1;
            return 0;
        }

        public int CompareTo(double v)
        {
            double t = h + l;
            if (t > v) return 1;
            if (t < v) return -1;
            v = v - h; t = l;
            if (t > v) return 1;
            if (t < v) return -1;
            return 0;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad operator +(Quad q, double v)
        {
            Quad r;
            r.h = (v + q.l) + q.h;
            r.l = q.h - r.h;
            v = (v + r.l) + q.l;
            r.l = q.h - (r.h + r.l) + v;
            return r;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad operator -(Quad q, double v)
        {
            Quad r;
            r.h = (q.l - v) + q.h;
            r.l = q.h - r.h;
            v = (r.l - v) + q.l;
            r.l = q.h - (r.h + r.l) + v;
            return r;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad operator +(Quad q, Quad v)
        {
            Quad r;
            double t = q.h + v.l;
            r.h = (q.l + v.h) + t;
            t = t - r.h;
            r.l = (t + v.h) + q.l;
            t = q.h - (r.h + t);
            r.l = (r.l + t) + v.l;
            return r;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad operator -(Quad q, Quad v)
        {
            Quad r;
            double t = q.h - v.l;
            r.h = (q.l - v.h) + t;
            t = t - r.h;
            r.l = (t - v.h) + q.l;
            t = q.h - (r.h + t);
            r.l = (r.l + t) - v.l;
            return r;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad Split(double v)
        {
            Quad r;
            r.h = (v - v * (1.0 * 0x7ffffff / 0x8000000)) * 0x8000000;
            r.l = v - r.h;
            return r;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad Split2(double v)
        {
            Quad r;
            r.l = v * 0x8000001;
            r.h = (v - r.l) + r.l;
            r.l = v - r.h;
            return r;
        }

        public static Quad Prod(double v, double w)
        {
            Quad r, vs = Split(v), ws = Split(w);
            r.h = v * w;
            r.l = vs.h * ws.h - r.h;
            r.l = vs.h * ws.l + r.l;
            r.l = vs.l * ws.h + r.l;
            r.l = vs.l * ws.l + r.l;
            return r;
        }

        public static Quad Sqr(double v)
        {
            Quad r, vs = Split(v);
            r.h = v * v;
            r.l = vs.h * vs.h - r.h;
            r.l = vs.h * vs.l + r.l;
            r.l = vs.l * vs.h + r.l;
            r.l = vs.l * vs.l + r.l;
            return r;
        }

        public static Quad Sqr(Quad q)
        {
            Quad r = Sqr(q.h);
            r.l += q.l * q.l;
            double t = q.h * q.l;
            r.FastAdd(t + t);
            return r;
        }

        public static Quad Sqrt(Quad q)
        {
            q.Scale(0.25);
            double t = Math.Sqrt(q.h);
            //Quad r = Sqr(t);
            //r.l = (q.h - r.h - r.l + q.l) / (t + t);
            //r.h = t;
            return q / t + t;
        }

        public static Quad operator ~(Quad q)
        {
            double t = 1 / q.h;
            Quad p = Prod(q.h, t);
            p.h = p.h - 1;
            p.FastAdd(q.l / q.h);
            return (p * -t) + t;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quad operator -(Quad q)
        {
            q.h = -q.h;
            q.l = -q.l;
            return q;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void FastAdd(double v)
        {
            double t = (l + v) + h;
            l += (h - t) + v;
            h = t;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Scale(double v)
        {
            h *= v;
            l *= v;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double Sign()
        {
            double r = this.h;
            if (r > 0) return 1;
            if (r < 0) return -1;
            return r;
        }

        public static Quad Abs(Quad x)
        {
            if (x.h < 0.0)
            {
                return -x;
            }
            return x;
        }

        public void FMA(double v, ref double w)
        {
            Quad t = Split(this.l);
            t.Scale(v);
            t.FastAdd(w);
            Quad r = Split(this.h);
            r.Scale(v);
            double s = (t.h + r.l) + r.h;
            this.h = s;
            t.FastAdd((r.h - s) + r.l);
            this.l = t.h;
            w = t.l;
        }

        public static Quad operator *(Quad q, double v)
        {
            Quad r = Prod(q.h, v);
            r.FastAdd(q.l * v);
            return r;
        }

        public static Quad operator /(Quad q, double v)
        {
            double t = q.h / v;
            Quad r = Prod(t, v);
            double s = (q.h - r.h) + q.l - r.l;
            r.h = t;
            r.l = s / v;
            return r;
        }

        public static Quad operator |(Quad q, Quad v)
        {
            Quad r = Prod(q.h, v.h);
            r.FastAdd(q.h * v.l);
            r.FastAdd(q.l * v.h);
            r.FastAdd(q.l * v.l);
            return r;
        }

        public static Quad operator *(Quad q, Quad v)
        {
            double t, s, p = q.l * v.l;
            Quad r = Prod(q.h, v.h);
            Quad r1 = Prod(q.h, v.l);
            s = (p + r1.l) + r.l;
            t = (r.l - s) + r1.l;
            r.l = s;
            p = p + t;
            s = (r.l + r1.h) + r.h;
            t = (r.h - s) + r1.h;
            r.h = s;
            r.l = r.l + t;
            Quad r2 = Prod(q.l, v.h);
            s = (p + r2.l) + r.l;
            t = (r.l - s) + r2.l;
            r.l = s;
            p = p + t;
            s = (r.l + r2.h) + r.h;
            t = (r.h - s) + r2.h;
            r.h = s;
            r.l = (r.l + t) + p;
            return r;
        }

        public static Quad operator &(Quad q, Quad v)
        {
            Quad r;
            double t = q.l, s = v.l;
            r.h = q.h * v.h;
            q = Split(q.h);
            q.l = q.l + t;
            v = Split(v.h);
            v.l = v.l + s;
            r.l = q.h * v.h - r.h;

            t = q.h * v.l + q.l * v.h;
            s = (r.l + t) + r.h;
            t = (r.h - s) + t + r.l;

            r.l = q.l * v.l;
            r.h = (r.l + t) + s;
            r.l = (s - r.h) + t + r.l;
            return r;
        }

        public static Quad Pow2(Quad q)
        {
            double v = Math.Round(q.h + q.l);
            if (v > 1024) return PositiveInfinity;
            else if (v <= -1074)
            {
                return q.h <= -1075 ? Zero : Epsilon;
            }

            q = PowHelp(new Quad(q.h - v, q.l));

            if (v == 1024) // edge case for values near MaxValue
            {
                if (q.h >= 1) return PositiveInfinity;
                q.Scale(2);
                v--;
            }
            q.Scale(Pow2((int)v));
            return q;
        }

        public static Quad Log2(Quad q)
        {
            if (q.h <= 0) return NegativeInfinity;

            q.Norm(out int n);

            return new Quad(n, LogHelp2(q));
        }

        public static bool IsFinite(double v)
        {
            return (BitConverter.DoubleToInt64Bits(v) & 0x7ff0000000000000) < 0x7ff0000000000000;
        }

        public static double Log2(double v)
        {
            // 1.66878860204731E+223 error -0.0860713320561217
            long t = BitConverter.DoubleToInt64Bits(v) & long.MaxValue;
            /*long m = t & 0xfffffffffffff;t <<= 1;
            long n = m >> 21;
            n *= n;
            n = n >> 11;
            if (m + n < 0x8000000000000)
                t = t + (n<<1);
            else
                t = t + n;
            t -= 0x7fe0000000000000;*/
            return (t - 0x3ff0000000000000) * 2.2204460492503131e-16;
        }

        public static double Pow2(int exp)
        {
            if (exp >= 1024)
            {
                return double.PositiveInfinity;
            }
            if (exp < -1074)
            {
                return 0;
            }
            long t;
            exp += 1023;
            if (exp > 0)
            {
                t = (long)exp << 52;
            }
            else
            {
                t = 0x8000000000000 >> -exp;
            }
            return BitConverter.Int64BitsToDouble(t);
        }

        public void Norm(out int exp)
        {
            long t = BitConverter.DoubleToInt64Bits(h);
            exp = (int)(t >> 52) & 2047;
            long m = t & 0xfffffffffffff;
            if (exp == 0)
            {
                if (m == 0)
                {
                    return;
                }
                while (true)
                {
                    t += m;
                    if (m >= 0x8000000000000)
                        break;
                    m += m;
                    exp--;
                }
                m = t & 0xfffffffffffff;
            }

            t |= 0x3ff0000000000000;
            if (m >= 0x6a09e667f3bcd)
            {
                m = ~0x4010000000000000;
                exp++;
            }
            else
            {
                m = ~0x4000000000000000;
            }
            h = BitConverter.Int64BitsToDouble(t & m);
            if (exp <= 0 || l == 0)
            {
                exp -= 1023;
                return;
            }

            exp = 2046 - exp;
            if (exp > 0)
            {
                t = (long)exp << 52;
            }
            else
            {
                t = 0x8000000000000 >> -exp;
            }
            exp = 1023 - exp;
            l *= BitConverter.Int64BitsToDouble(t);
        }

        public static Quad Sin(Quad q)
        {
            return SinCos(q, 0);
        }

        public static Quad Cos(Quad q)
        {
            return SinCos(q, 1);
        }

        // maximal input range -PI*2^102 - PI*2^102
        public static Quad SinCos(Quad q, int quadrant)
        {
            Quad d = q * ~PI;
            double h = d.h + 1.5 * (1L << 52) * (1L << 51) - 1.5 * (1L << 52) * (1L << 51);
            double l = (d.h - h) + 1.5 * (1L << 51) + d.l;
            long t = BitConverter.DoubleToInt64Bits(l) + quadrant;
            l -= 1.5 * (1L << 51);

            q = q - Prod(h, PI.h) - Prod(l, PI.h);
            q = q - Prod(h, PI.l) - Prod(l, PI.l);
            q = q - Prod(h, pilow.h) - Prod(l, pilow.h);
            q = q - Prod(h, pilow.l) - Prod(l, pilow.l);
            q = ((t & 1) == 0) ? SinHelp(q) : CosHelp(q);
            return ((t & 2) == 0) ? q : -q;
        }

        // slow convergence
        // optimal input range -1/3 - 1/3
        private static Quad LogHelp1(Quad q)
        {
            Quad q1 = -q;
            Quad p = q *= inv2ln2;
            for (int i = 2; i <= 67; i++)
            {
                p *= q1;
                Quad t = q + p / i;
                if (t.h == q.h && t.l == q.l)
                    break;
                q = t;
            }
            return q;
        }

        // fast convergence
        // optimal input range : 1/sqrt(2) - sqrt(2)
        private static Quad LogHelp2(Quad q)
        {
            double s = q.Sign();
            q = new Quad(q.h - s, q.l) * ~(q + s);

            Quad q2 = Sqr(q);
            Quad p = q *= inv2ln2;
            for (int i = 3; i <= 43; i+=2)
            {
                p = p * q2;
                q += p / i;
            }
            return q;
        }

        // optimal input range : -tan(pi/16) - tan(pi/16)
        private static Quad AtanHelp(Quad q)
        {
            Quad q2 = -Sqr(q);
            Quad p = q;
            for (int i = 3; i <= 43; i+=2)
            {
                p = p * q2;
                q += p / i;
            }
            return q;
        }

        private static Quad SinHelp(Quad q)
        {
            Quad q2 = Sqr(q);
            q2.Scale(-0.5);
            Quad p = q;
            for (int i = 1; i <= 13; i++)
            {
                p = p * q2 / (i * ((i << 1) + 1));
                q += p;
            }
            return q;
        }

        private static Quad CosHelp(Quad q)
        {
            Quad q2 = Sqr(q);
            q2.Scale(-0.5);
            Quad p = q = q2;
            for (int i = 2; i <= 14; i++)
            {
                p = p * q2 / (i * ((i << 1) - 1));
                q += p;
            }
            return new Quad(1, q);
        }

        // optimal input range : -0.5 - +0.5
        private static Quad PowHelp(Quad q)
        {
            Quad p = q *= ln2;
            Quad q1 = q;
            for (int i = 2; i <= 23; i++)
            {
                p = p * q1 / i;
                q += p;
            }
            return new Quad(1, q);
        }

        public static Quad Pow(double x, int n)
        {
            Quad q = x;
            if (n < 0)
            {
                return ~Pow(q, -n);
            }
            return Pow(q, n);
        }

        public static Quad Pow(Quad q, int n)
        {
            Quad r;
            if (n < 0)
            {
                q = ~q;
                n = -n;
            }
            if ((n & 1) != 0)
            {
                r = q;
            }
            else
            {
                r = q.Sign();
            }

            while ((n >>= 1) != 0)
            {
                q = Sqr(q);
                if ((n & 1) != 0)
                {
                    r *= q;
                }
            }
            return r;
        }

        public static Quad Parse(string s)
        {
            return Parse(s, NumberFormatInfo.CurrentInfo);
        }

        public static Quad Parse(string s, IFormatProvider provider)
        {
            return Parse(s, NumberFormatInfo.GetInstance(provider));
        }

        public static Quad Parse(string s, NumberFormatInfo info)
        {
            if (s == null)
            {
                throw new ArgumentNullException("s");
            }
            if (!TryParse(s, info, out Quad r))
            {
                if (IsFinite(r.h))
                {
                    throw new FormatException();
                }
                throw new ArgumentOutOfRangeException("s");
            }
            return r;
        }

        public static bool TryParse(string s, out Quad q)
        {
            return TryParse(s, NumberFormatInfo.CurrentInfo, out q);
        }

        public static bool TryParse(string s, IFormatProvider provider, out Quad q)
        {
            return TryParse(s, NumberFormatInfo.GetInstance(provider), out q);
        }

        public static bool TryParse(string s, NumberFormatInfo info, out Quad q)
        {
            q = new Quad();
            if (s == null)
            {
                return false;
            }
            int end = s.Length;
            do
            {
                if (end == 0)
                {
                    return false;
                }
            }
            while (char.IsWhiteSpace(s[--end]));
            int start = 0;
            int exp, tmp, e;
            char sign;
            while (char.IsWhiteSpace(sign = s[start]))
            {
                start++;
            }
            end++;

            string n = info.PositiveInfinitySymbol;
            tmp = ',' - sign;
            if (tmp == 1 || tmp == -1)
            {
                if (++start == end)
                {
                    return false;
                }
            }
            else if (sign == n[0])
            {
                tmp = 1;
            }
            else
            {
                tmp = 0;
                n = info.NaNSymbol;
            }

            if (n.Length == end - start && string.CompareOrdinal(n, 0, s, start, end - start) == 0)
            {
                q.h = double.PositiveInfinity * tmp;
                return true;
            }

            exp = 1;
            double w = 0, v = 0;
            do
            {
                tmp = s[start] - '0';
                if (tmp < 0)
                {
                    if (exp == 0)
                    {
                        exp = '-' - '0' - tmp;
                        if (exp == -1) // point
                        {
                            continue;
                        }
                        if (exp == 1) // comma
                        {
                            continue;
                        }
                    }
                    return false;
                }
                else if (tmp > 9)
                {
                    if (tmp != 'E' - '0' && tmp != 'e' - '0')
                    {
                        return false;
                    }
                    if (++start == end)
                    {
                        return false;
                    }
                    break;
                }
                else if (exp != 0)
                {
                    exp--;
                }
                if (w == 0)
                {
                    if (tmp == 0)
                    {
                        continue;
                    }
                    w = 1;
                }
                v = v * 10 + tmp;
                if ((w *= 10) >= 1E11)
                {
                    q.FMA(w, ref v);
                    w = 1;
                }
            }
            while (++start < end);
            q.FMA(w, ref v);

            if (exp >= 0)
            {
                exp = -1;
            }
            if (start < end) // has exponent
            {
                tmp = ',' - s[start];
                if (tmp == 1 || tmp == -1)
                {
                    exp *= tmp;
                    if (++start == end)
                    {
                        return false;
                    }
                }

                e = exp;
                exp = 0;
                do
                {
                    tmp = s[start] - '0';
                    if (tmp < 0 || tmp > 9)
                    {
                        return false;
                    }
                    // prevent int overflow
                    if (exp <= 107374183)
                    {
                        exp = exp * 10 + tmp;
                    }
                }
                while (++start < end);
                exp += e;
                if (e >= 0)
                {
                    exp = -exp;
                }
            }
            if (q.h != 0)
            {
                q *= Pow(sign == '-' ? -5 : 5, exp + 1);
                q.Scale(Pow2(exp + 1));
            }
            return IsFinite(q.h);
        }

        public string ToStringInvariant()
        {
            return ToString(NumberFormatInfo.InvariantInfo);
        }

        public override string ToString()
        {
            return ToString(NumberFormatInfo.CurrentInfo);
        }

        public string ToString(IFormatProvider provider)
        {
            return ToString(NumberFormatInfo.GetInstance(provider));
        }

        public string ToString(NumberFormatInfo info)
        {
            if (this.h == 0.0)
            {
                return this.l.ToString(info);
            }
            if (!IsFinite(this.h))
            {
                return this.h > 0 ? info.PositiveInfinitySymbol : this.h < 0 ? info.NegativeInfinitySymbol : info.NaNSymbol;
            }

            double v = Log2(this.h);
            if (v > 0)
                v = v + 0.0860713320561217;
            int exp = (int)(v * 0.3010299956639812);
            int len = 0;
            string s = null;
            Quad q = this;
            q.Scale(Pow2(10 - exp));
            q *= Pow(this.h < 0 ? -5 : 5, 10 - exp);

            v = .5e-24;
            // round the last digit
            q.FMA(1, ref v);
            // check for underestimate
            if (q.CompareTo(1E10) < 0)
            {
                v *= 10;
                q.FMA(10, ref v);
                exp--;
            }
            if (exp > -5 && exp < 0)
            {
                s = new string('0', -exp);
                exp = 0;
            }

            do
            {
                long t = q.ToLong();
                q = new Quad(q.h - t, q.l);
                v *= 1E11;
                q.FMA(1E11, ref v);
                s += t.ToString("00000000000", info);
            }
            while (++len < 3);
            len = s.Length;
            do if (--len == exp && exp < 16) break;
            while (s[len] == '0');
            s = s.Substring(0, len + 1);

            if (exp >= 16 || exp < 0)
            {
                s += exp.ToString("\\E+00;\\E-00", info);
                exp = 0;
            }
            if (exp < len)
            {
                s = s.Insert(exp + 1, info.NumberDecimalSeparator);
            }
            if (this.h < 0.0)
            {
                s = info.NegativeSign + s;
            }
            return s;
        }
    }
}
