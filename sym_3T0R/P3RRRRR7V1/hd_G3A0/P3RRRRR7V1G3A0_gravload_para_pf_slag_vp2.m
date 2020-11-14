% Calculate Gravitation load for parallel robot
% P3RRRRR7V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:53:18
% EndTime: 2020-08-07 08:53:20
% DurationCPUTime: 1.46s
% Computational Cost: add. (855->236), mult. (1287->320), div. (60->14), fcn. (951->60), ass. (0->148)
t1827 = sin(qJ(1,3));
t1844 = -pkin(5) - pkin(4);
t1795 = t1827 * t1844;
t1836 = cos(qJ(1,3));
t1834 = cos(qJ(3,3));
t1918 = t1834 * pkin(2);
t1790 = pkin(1) + t1918;
t1835 = cos(qJ(2,3));
t1825 = sin(qJ(3,3));
t1826 = sin(qJ(2,3));
t1887 = t1825 * t1826;
t1862 = pkin(2) * t1887 - t1790 * t1835;
t1927 = t1836 * t1862 + t1795;
t1830 = sin(qJ(1,2));
t1796 = t1830 * t1844;
t1839 = cos(qJ(1,2));
t1837 = cos(qJ(3,2));
t1917 = t1837 * pkin(2);
t1792 = pkin(1) + t1917;
t1838 = cos(qJ(2,2));
t1828 = sin(qJ(3,2));
t1829 = sin(qJ(2,2));
t1885 = t1828 * t1829;
t1861 = pkin(2) * t1885 - t1792 * t1838;
t1926 = t1839 * t1861 + t1796;
t1833 = sin(qJ(1,1));
t1797 = t1833 * t1844;
t1842 = cos(qJ(1,1));
t1840 = cos(qJ(3,1));
t1916 = t1840 * pkin(2);
t1794 = pkin(1) + t1916;
t1841 = cos(qJ(2,1));
t1831 = sin(qJ(3,1));
t1832 = sin(qJ(2,1));
t1883 = t1831 * t1832;
t1860 = pkin(2) * t1883 - t1794 * t1841;
t1925 = t1842 * t1860 + t1797;
t1915 = m(3) * pkin(1) + mrSges(2,1);
t1771 = mrSges(3,1) * t1840 - mrSges(3,2) * t1831 + t1915;
t1783 = mrSges(3,1) * t1831 + mrSges(3,2) * t1840 + mrSges(2,2);
t1924 = t1841 * t1771 - t1783 * t1832;
t1770 = mrSges(3,1) * t1837 - mrSges(3,2) * t1828 + t1915;
t1782 = mrSges(3,1) * t1828 + mrSges(3,2) * t1837 + mrSges(2,2);
t1923 = t1838 * t1770 - t1782 * t1829;
t1769 = mrSges(3,1) * t1834 - mrSges(3,2) * t1825 + t1915;
t1781 = mrSges(3,1) * t1825 + mrSges(3,2) * t1834 + mrSges(2,2);
t1922 = t1835 * t1769 - t1781 * t1826;
t1921 = 0.2e1 * pkin(2);
t1920 = 0.2e1 * t1844;
t1914 = -qJ(3,1) + qJ(1,1);
t1913 = qJ(3,1) + qJ(1,1);
t1912 = -qJ(3,2) + qJ(1,2);
t1911 = qJ(3,2) + qJ(1,2);
t1910 = -qJ(3,3) + qJ(1,3);
t1909 = qJ(3,3) + qJ(1,3);
t1908 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1907 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1906 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1905 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1904 = 0.2e1 * qJ(2,3) + qJ(1,3);
t1903 = -0.2e1 * qJ(2,3) + qJ(1,3);
t1822 = legFrame(3,2);
t1804 = sin(t1822);
t1807 = cos(t1822);
t1772 = g(1) * t1804 + g(2) * t1807;
t1819 = qJ(2,3) + qJ(3,3);
t1798 = sin(t1819);
t1801 = cos(t1819);
t1810 = 0.1e1 / t1825;
t1775 = g(1) * t1807 - g(2) * t1804;
t1865 = g(3) * t1827 - t1775 * t1836;
t1902 = ((-mrSges(3,1) * t1772 - mrSges(3,2) * t1865) * t1801 + t1798 * (-mrSges(3,1) * t1865 + mrSges(3,2) * t1772)) * t1810;
t1823 = legFrame(2,2);
t1805 = sin(t1823);
t1808 = cos(t1823);
t1773 = g(1) * t1805 + g(2) * t1808;
t1820 = qJ(2,2) + qJ(3,2);
t1799 = sin(t1820);
t1802 = cos(t1820);
t1811 = 0.1e1 / t1828;
t1776 = g(1) * t1808 - g(2) * t1805;
t1864 = g(3) * t1830 - t1776 * t1839;
t1901 = ((-mrSges(3,1) * t1773 - mrSges(3,2) * t1864) * t1802 + t1799 * (-mrSges(3,1) * t1864 + mrSges(3,2) * t1773)) * t1811;
t1824 = legFrame(1,2);
t1806 = sin(t1824);
t1809 = cos(t1824);
t1774 = g(1) * t1806 + g(2) * t1809;
t1821 = qJ(2,1) + qJ(3,1);
t1800 = sin(t1821);
t1803 = cos(t1821);
t1812 = 0.1e1 / t1831;
t1777 = g(1) * t1809 - g(2) * t1806;
t1863 = g(3) * t1833 - t1777 * t1842;
t1900 = ((-mrSges(3,1) * t1774 - mrSges(3,2) * t1863) * t1803 + t1800 * (-mrSges(3,1) * t1863 + mrSges(3,2) * t1774)) * t1812;
t1788 = m(3) * pkin(4) - mrSges(1,2) + mrSges(2,3) + mrSges(3,3);
t1787 = g(3) * t1788;
t1843 = g(3) * mrSges(1,1);
t1899 = 0.1e1 / (pkin(1) * t1835 + pkin(2) * t1801) * ((g(3) * t1922 + t1843) * t1836 + t1827 * t1787 + (-t1788 * t1836 + t1827 * (mrSges(1,1) + t1922)) * t1775);
t1898 = 0.1e1 / (pkin(1) * t1838 + pkin(2) * t1802) * ((g(3) * t1923 + t1843) * t1839 + t1830 * t1787 + (-t1788 * t1839 + t1830 * (mrSges(1,1) + t1923)) * t1776);
t1897 = 0.1e1 / (pkin(1) * t1841 + pkin(2) * t1803) * ((g(3) * t1924 + t1843) * t1842 + t1833 * t1787 + (-t1788 * t1842 + t1833 * (mrSges(1,1) + t1924)) * t1777);
t1813 = t1834 ^ 2;
t1784 = pkin(1) * t1834 + t1813 * t1921 - pkin(2);
t1893 = t1784 * t1836;
t1815 = t1837 ^ 2;
t1785 = pkin(1) * t1837 + t1815 * t1921 - pkin(2);
t1892 = t1785 * t1839;
t1817 = t1840 ^ 2;
t1786 = pkin(1) * t1840 + t1817 * t1921 - pkin(2);
t1891 = t1786 * t1842;
t1890 = (pkin(1) + 0.2e1 * t1918) * t1825;
t1889 = (pkin(1) + 0.2e1 * t1917) * t1828;
t1888 = (pkin(1) + 0.2e1 * t1916) * t1831;
t1886 = t1826 * t1784;
t1884 = t1829 * t1785;
t1882 = t1832 * t1786;
t1877 = t1825 * t1918;
t1876 = t1828 * t1917;
t1875 = t1831 * t1916;
t1745 = (-t1769 * t1865 + t1772 * t1781) * t1826 - t1835 * (t1772 * t1769 + t1781 * t1865);
t1874 = t1745 / t1862 * t1810;
t1746 = (-t1770 * t1864 + t1773 * t1782) * t1829 - t1838 * (t1770 * t1773 + t1782 * t1864);
t1873 = t1746 / t1861 * t1811;
t1747 = (-t1771 * t1863 + t1774 * t1783) * t1832 - t1841 * (t1774 * t1771 + t1783 * t1863);
t1872 = t1747 / t1860 * t1812;
t1871 = t1836 * t1887;
t1870 = t1827 * t1899;
t1869 = t1839 * t1885;
t1868 = t1830 * t1898;
t1867 = t1842 * t1883;
t1866 = t1833 * t1897;
t1859 = pkin(1) * t1871 + (t1871 * t1921 + t1795) * t1834;
t1858 = pkin(1) * t1869 + (t1869 * t1921 + t1796) * t1837;
t1857 = pkin(1) * t1867 + (t1867 * t1921 + t1797) * t1840;
t1856 = 0.1e1 / pkin(1);
t1855 = 0.1e1 / pkin(2);
t1852 = 0.2e1 * qJ(3,1);
t1849 = 0.2e1 * qJ(3,2);
t1846 = 0.2e1 * qJ(3,3);
t1818 = t1841 ^ 2;
t1816 = t1838 ^ 2;
t1814 = t1835 ^ 2;
t1762 = pkin(2) * t1831 * t1841 + t1794 * t1832;
t1761 = pkin(2) * t1828 * t1838 + t1792 * t1829;
t1760 = pkin(2) * t1825 * t1835 + t1790 * t1826;
t1756 = -t1797 * t1883 + (t1817 - 0.1e1) * t1842 * pkin(2);
t1755 = -t1796 * t1885 + (t1815 - 0.1e1) * t1839 * pkin(2);
t1754 = -t1795 * t1887 + (t1813 - 0.1e1) * t1836 * pkin(2);
t1 = [-t1807 * t1870 - t1808 * t1868 - t1809 * t1866 - g(1) * m(4) + (-((t1806 * t1888 + t1809 * t1891) * t1818 + (t1806 * t1882 - t1809 * t1857) * t1841 - t1756 * t1809 - t1806 * t1875) * t1872 - ((t1805 * t1889 + t1808 * t1892) * t1816 + (t1805 * t1884 - t1808 * t1858) * t1838 - t1755 * t1808 - t1805 * t1876) * t1873 - ((t1804 * t1890 + t1807 * t1893) * t1814 + (t1804 * t1886 - t1807 * t1859) * t1835 - t1754 * t1807 - t1804 * t1877) * t1874 + ((-t1806 * t1762 + t1809 * t1925) * t1900 + (-t1805 * t1761 + t1808 * t1926) * t1901 + (-t1804 * t1760 + t1807 * t1927) * t1902) * t1855) * t1856; t1804 * t1870 + t1805 * t1868 + t1806 * t1866 - g(2) * m(4) + (-((-t1806 * t1891 + t1809 * t1888) * t1818 + (t1806 * t1857 + t1809 * t1882) * t1841 + t1756 * t1806 - t1809 * t1875) * t1872 - ((-t1805 * t1892 + t1808 * t1889) * t1816 + (t1805 * t1858 + t1808 * t1884) * t1838 + t1755 * t1805 - t1808 * t1876) * t1873 - ((-t1804 * t1893 + t1807 * t1890) * t1814 + (t1804 * t1859 + t1807 * t1886) * t1835 + t1754 * t1804 - t1807 * t1877) * t1874 + ((-t1809 * t1762 - t1806 * t1925) * t1900 + (-t1808 * t1761 - t1805 * t1926) * t1901 + (-t1807 * t1760 - t1804 * t1927) * t1902) * t1855) * t1856; -g(3) * m(4) - t1836 * t1899 - t1839 * t1898 - t1842 * t1897 + (((cos(-qJ(2,1) + t1914) + cos(qJ(2,1) + t1913)) * t1920 + (sin(-0.2e1 * qJ(3,1) + t1907) + sin(t1852 + t1908) + 0.2e1 * t1833) * pkin(2) + (sin(-qJ(3,1) + t1907) + sin(qJ(3,1) + t1908) + sin(t1914) + sin(t1913)) * pkin(1)) / ((-sin(t1852 + qJ(2,1)) + t1832) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t1800) * pkin(1)) * t1747 + ((cos(-qJ(2,2) + t1912) + cos(qJ(2,2) + t1911)) * t1920 + (sin(-0.2e1 * qJ(3,2) + t1905) + sin(t1849 + t1906) + 0.2e1 * t1830) * pkin(2) + (sin(-qJ(3,2) + t1905) + sin(qJ(3,2) + t1906) + sin(t1912) + sin(t1911)) * pkin(1)) / ((-sin(t1849 + qJ(2,2)) + t1829) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t1799) * pkin(1)) * t1746 + ((cos(-qJ(2,3) + t1910) + cos(qJ(2,3) + t1909)) * t1920 + (sin(-0.2e1 * qJ(3,3) + t1903) + sin(t1846 + t1904) + 0.2e1 * t1827) * pkin(2) + (sin(-qJ(3,3) + t1903) + sin(qJ(3,3) + t1904) + sin(t1910) + sin(t1909)) * pkin(1)) / ((-sin(t1846 + qJ(2,3)) + t1826) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t1798) * pkin(1)) * t1745) * t1856 / 0.2e1 + ((-t1833 * t1860 + t1842 * t1844) * t1900 + (-t1830 * t1861 + t1839 * t1844) * t1901 + (-t1827 * t1862 + t1836 * t1844) * t1902) * t1855 * t1856;];
taugX  = t1;
