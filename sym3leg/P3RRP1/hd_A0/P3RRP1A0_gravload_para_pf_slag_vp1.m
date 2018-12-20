% Calculate Gravitation load for parallel robot
% P3RRP1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m [4x1]
%   mass of all robot links (including platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 18:08
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3RRP1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 18:00:28
% EndTime: 2018-12-20 18:00:30
% DurationCPUTime: 2.57s
% Computational Cost: add. (2328->370), mult. (4797->556), div. (36->3), fcn. (2138->20), ass. (0->252)
t1884 = (qJ(3,3) ^ 2);
t1894 = pkin(1) ^ 2;
t1876 = cos(qJ(1,3));
t1857 = t1876 ^ 2;
t1928 = t1857 * t1894;
t1981 = (2 * t1884) + t1894 - t1928;
t1885 = (qJ(3,2) ^ 2);
t1877 = cos(qJ(1,2));
t1858 = t1877 ^ 2;
t1927 = t1858 * t1894;
t1980 = (2 * t1885) + t1894 - t1927;
t1886 = (qJ(3,1) ^ 2);
t1878 = cos(qJ(1,1));
t1859 = t1878 ^ 2;
t1926 = t1859 * t1894;
t1979 = (2 * t1886) + t1894 - t1926;
t1978 = 2 * pkin(1);
t1870 = legFrame(3,3);
t1840 = sin(t1870);
t1822 = t1840 * qJ(3,3);
t1819 = pkin(2) * t1822;
t1977 = -0.2e1 * t1819;
t1871 = legFrame(2,3);
t1841 = sin(t1871);
t1823 = t1841 * qJ(3,2);
t1820 = pkin(2) * t1823;
t1976 = -0.2e1 * t1820;
t1872 = legFrame(1,3);
t1842 = sin(t1872);
t1824 = t1842 * qJ(3,1);
t1821 = pkin(2) * t1824;
t1975 = -0.2e1 * t1821;
t1873 = sin(qJ(1,3));
t1974 = -0.2e1 * t1873;
t1874 = sin(qJ(1,2));
t1973 = -0.2e1 * t1874;
t1875 = sin(qJ(1,1));
t1972 = -0.2e1 * t1875;
t1971 = -0.2e1 * t1876;
t1970 = -0.2e1 * t1877;
t1969 = -0.2e1 * t1878;
t1968 = m(1) * rSges(1,2);
t1967 = pkin(2) * qJ(3,1);
t1966 = pkin(2) * qJ(3,2);
t1965 = pkin(2) * qJ(3,3);
t1843 = cos(t1870);
t1964 = pkin(2) * t1843;
t1844 = cos(t1871);
t1963 = pkin(2) * t1844;
t1845 = cos(t1872);
t1962 = pkin(2) * t1845;
t1961 = pkin(2) * t1876;
t1960 = pkin(2) * t1877;
t1959 = pkin(2) * t1878;
t1825 = t1840 * pkin(2);
t1826 = t1841 * pkin(2);
t1827 = t1842 * pkin(2);
t1958 = qJ(3,1) * t1845;
t1957 = qJ(3,1) * t1875;
t1956 = qJ(3,2) * t1844;
t1955 = qJ(3,2) * t1874;
t1954 = qJ(3,3) * t1843;
t1953 = qJ(3,3) * t1873;
t1860 = qJ(1,3) + qJ(2,3);
t1834 = sin(t1860);
t1952 = t1834 * qJ(3,3);
t1861 = qJ(1,2) + qJ(2,2);
t1835 = sin(t1861);
t1951 = t1835 * qJ(3,2);
t1862 = qJ(1,1) + qJ(2,1);
t1836 = sin(t1862);
t1950 = t1836 * qJ(3,1);
t1837 = cos(t1860);
t1828 = t1837 ^ 2;
t1893 = pkin(2) ^ 2;
t1846 = -t1884 + t1893;
t1831 = 0.1e1 + t1846;
t1863 = 0.2e1 + t1894;
t1925 = t1873 * t1876;
t1895 = -t1831 * t1857 + 0.2e1 * t1925 * t1965;
t1922 = -0.1e1 / 0.2e1 - t1893 / 0.2e1;
t1748 = 0.1e1 / (((t1831 * t1925 + (0.2e1 * t1857 - 0.1e1) * t1965) * pkin(1) * t1834 + t1953) * t1837 * t1978 + pkin(1) * t1952 * t1971 - t1884 * t1863 + (0.2e1 * (t1884 / 0.2e1 - t1895 + t1922) * t1828 + t1895) * t1894);
t1796 = -g(1) * t1840 + g(2) * t1843;
t1799 = g(1) * t1843 + g(2) * t1840;
t1760 = -t1796 * t1837 + t1799 * t1834;
t1949 = t1748 * t1760;
t1838 = cos(t1861);
t1829 = t1838 ^ 2;
t1847 = -t1885 + t1893;
t1832 = 0.1e1 + t1847;
t1924 = t1874 * t1877;
t1896 = -t1832 * t1858 + 0.2e1 * t1924 * t1966;
t1749 = 0.1e1 / (((t1832 * t1924 + (0.2e1 * t1858 - 0.1e1) * t1966) * pkin(1) * t1835 + t1955) * t1838 * t1978 + pkin(1) * t1951 * t1970 - t1885 * t1863 + (0.2e1 * (t1885 / 0.2e1 - t1896 + t1922) * t1829 + t1896) * t1894);
t1797 = -g(1) * t1841 + g(2) * t1844;
t1800 = g(1) * t1844 + g(2) * t1841;
t1761 = -t1797 * t1838 + t1800 * t1835;
t1948 = t1749 * t1761;
t1839 = cos(t1862);
t1830 = t1839 ^ 2;
t1848 = -t1886 + t1893;
t1833 = 0.1e1 + t1848;
t1923 = t1875 * t1878;
t1897 = -t1833 * t1859 + 0.2e1 * t1923 * t1967;
t1750 = 0.1e1 / (((t1833 * t1923 + (0.2e1 * t1859 - 0.1e1) * t1967) * pkin(1) * t1836 + t1957) * t1839 * t1978 + pkin(1) * t1950 * t1969 - t1886 * t1863 + (0.2e1 * (t1886 / 0.2e1 - t1897 + t1922) * t1830 + t1897) * t1894);
t1798 = -g(1) * t1842 + g(2) * t1845;
t1801 = g(1) * t1845 + g(2) * t1842;
t1762 = -t1798 * t1839 + t1801 * t1836;
t1947 = t1750 * t1762;
t1784 = t1822 + t1964;
t1946 = t1784 * t1873;
t1788 = t1823 + t1963;
t1945 = t1788 * t1874;
t1792 = t1824 + t1962;
t1944 = t1792 * t1875;
t1943 = t1822 * t1837;
t1942 = t1823 * t1838;
t1941 = t1824 * t1839;
t1940 = t1831 * t1843;
t1939 = t1832 * t1844;
t1938 = t1833 * t1845;
t1937 = t1834 * t1837;
t1936 = t1834 * t1840;
t1935 = t1835 * t1838;
t1934 = t1836 * t1839;
t1933 = t1836 * t1842;
t1932 = t1837 * t1843;
t1931 = t1839 * t1845;
t1930 = t1841 * t1835;
t1929 = t1844 * t1838;
t1868 = t1893 / 0.2e1;
t1921 = t1868 - t1884 / 0.2e1;
t1920 = t1868 - t1885 / 0.2e1;
t1919 = t1868 - t1886 / 0.2e1;
t1918 = pkin(2) * t1958;
t1917 = pkin(2) * t1956;
t1916 = pkin(2) * t1954;
t1915 = t1843 * t1952;
t1914 = t1844 * t1951;
t1913 = t1845 * t1950;
t1912 = t1894 * t1925;
t1911 = t1894 * t1924;
t1910 = t1894 * t1923;
t1909 = -(2 * t1884) - t1928;
t1908 = -(2 * t1885) - t1927;
t1907 = -(2 * t1886) - t1926;
t1906 = t1843 * t1912;
t1905 = t1844 * t1911;
t1904 = t1845 * t1910;
t1785 = t1822 - t1964;
t1786 = t1825 + t1954;
t1903 = -t1785 * t1857 - t1786 * t1925;
t1902 = -t1785 * t1925 + t1786 * t1857;
t1789 = t1823 - t1963;
t1790 = t1826 + t1956;
t1901 = -t1789 * t1858 - t1790 * t1924;
t1900 = -t1789 * t1924 + t1790 * t1858;
t1793 = t1824 - t1962;
t1794 = t1827 + t1958;
t1899 = -t1793 * t1859 - t1794 * t1923;
t1898 = -t1793 * t1923 + t1794 * t1859;
t1892 = koppelP(1,1);
t1891 = koppelP(2,1);
t1890 = koppelP(3,1);
t1889 = koppelP(1,2);
t1888 = koppelP(2,2);
t1887 = koppelP(3,2);
t1883 = rSges(4,1);
t1882 = rSges(4,2);
t1881 = xP(3);
t1880 = m(2) * rSges(2,2);
t1879 = pkin(2) + rSges(3,1);
t1867 = qJ(3,1) + rSges(3,3);
t1866 = qJ(3,2) + rSges(3,3);
t1865 = qJ(3,3) + rSges(3,3);
t1864 = 0.1e1 + t1894;
t1853 = cos(t1881);
t1852 = sin(t1881);
t1818 = m(2) * rSges(2,1) + m(3) * t1879;
t1817 = -0.2e1 * t1918;
t1816 = -0.2e1 * t1917;
t1815 = -0.2e1 * t1916;
t1814 = 0.2e1 * t1821;
t1813 = 0.2e1 * t1820;
t1812 = 0.2e1 * t1819;
t1811 = -m(3) * t1867 + t1880;
t1810 = -m(3) * t1866 + t1880;
t1809 = -m(3) * t1865 + t1880;
t1808 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1807 = -t1957 + t1959;
t1806 = -t1955 + t1960;
t1805 = -t1953 + t1961;
t1804 = t1842 * t1910;
t1803 = t1841 * t1911;
t1802 = t1840 * t1912;
t1795 = t1827 - t1958;
t1791 = t1826 - t1956;
t1787 = t1825 - t1954;
t1783 = -t1852 * t1889 + t1853 * t1892;
t1782 = -t1852 * t1888 + t1853 * t1891;
t1781 = -t1852 * t1887 + t1853 * t1890;
t1780 = -t1852 * t1892 - t1853 * t1889;
t1779 = -t1852 * t1891 - t1853 * t1888;
t1778 = -t1852 * t1890 - t1853 * t1887;
t1777 = t1842 * t1848 + t1817;
t1776 = t1841 * t1847 + t1816;
t1775 = t1840 * t1846 + t1815;
t1774 = t1833 * t1842 + t1817;
t1773 = t1832 * t1841 + t1816;
t1772 = t1831 * t1840 + t1815;
t1771 = t1774 * t1878;
t1770 = t1773 * t1877;
t1769 = t1772 * t1876;
t1768 = t1795 * t1878 - t1944;
t1767 = t1792 * t1878 + t1795 * t1875;
t1766 = t1791 * t1877 - t1945;
t1765 = t1788 * t1877 + t1791 * t1874;
t1764 = t1787 * t1876 - t1946;
t1763 = t1784 * t1876 + t1787 * t1873;
t1759 = (t1845 * t1848 + t1814) * t1878 + t1777 * t1875;
t1758 = (t1844 * t1847 + t1813) * t1877 + t1776 * t1874;
t1757 = (t1843 * t1846 + t1812) * t1876 + t1775 * t1873;
t1756 = (t1814 + t1938) * t1878 + t1875 * t1774;
t1755 = (t1813 + t1939) * t1877 + t1874 * t1773;
t1754 = (t1812 + t1940) * t1876 + t1873 * t1772;
t1753 = t1777 * t1878 + (t1845 * t1919 + t1821) * t1972;
t1752 = t1776 * t1877 + (t1844 * t1920 + t1820) * t1973;
t1751 = t1775 * t1876 + (t1843 * t1921 + t1819) * t1974;
t1747 = ((-t1798 * t1879 - t1801 * t1867) * m(3) + m(2) * (-rSges(2,1) * t1798 + rSges(2,2) * t1801)) * t1839 + t1836 * ((-t1798 * t1867 + t1801 * t1879) * m(3) + m(2) * (rSges(2,1) * t1801 + rSges(2,2) * t1798));
t1746 = ((-t1797 * t1879 - t1800 * t1866) * m(3) + m(2) * (-rSges(2,1) * t1797 + rSges(2,2) * t1800)) * t1838 + t1835 * ((-t1797 * t1866 + t1800 * t1879) * m(3) + m(2) * (rSges(2,1) * t1800 + rSges(2,2) * t1797));
t1745 = ((-t1796 * t1879 - t1799 * t1865) * m(3) + m(2) * (-rSges(2,1) * t1796 + rSges(2,2) * t1799)) * t1837 + t1834 * ((-t1796 * t1865 + t1799 * t1879) * m(3) + m(2) * (rSges(2,1) * t1799 + rSges(2,2) * t1796));
t1744 = (-t1931 + t1933) * qJ(3,1) + (-t1771 * t1830 - t1756 * t1934 + (t1842 * t1893 + t1842 - t1918) * t1878 + (-(t1975 - t1938) * t1830 - qJ(3,1) * t1795) * t1875) * pkin(1);
t1743 = (-t1929 + t1930) * qJ(3,2) + (-t1770 * t1829 - t1755 * t1935 + (t1841 * t1893 + t1841 - t1917) * t1877 + (-(t1976 - t1939) * t1829 - qJ(3,2) * t1791) * t1874) * pkin(1);
t1742 = (-t1932 + t1936) * qJ(3,3) + (-t1769 * t1828 - t1754 * t1937 + (t1840 * t1893 + t1840 - t1916) * t1876 + (-(t1977 - t1940) * t1828 - qJ(3,3) * t1787) * t1873) * pkin(1);
t1741 = -t1913 - t1941 + (t1756 * t1830 - (t1771 + ((0.1e1 / 0.2e1 + t1919) * t1845 + t1821) * t1972) * t1934 - (t1845 * t1893 + t1821 + t1845) * t1878 + qJ(3,1) * t1944) * pkin(1);
t1740 = -t1914 - t1942 + (t1755 * t1829 - (t1770 + ((0.1e1 / 0.2e1 + t1920) * t1844 + t1820) * t1973) * t1935 - (t1844 * t1893 + t1820 + t1844) * t1877 + qJ(3,2) * t1945) * pkin(1);
t1739 = -t1915 - t1943 + (t1754 * t1828 - (t1769 + ((0.1e1 / 0.2e1 + t1921) * t1843 + t1819) * t1974) * t1937 - (t1843 * t1893 + t1819 + t1843) * t1876 + qJ(3,3) * t1946) * pkin(1);
t1738 = (-t1798 * t1818 + t1801 * t1811) * t1839 + (t1798 * t1811 + t1801 * t1818) * t1836 + (-t1798 * t1808 + t1801 * t1968) * t1878 + t1875 * (t1798 * t1968 + t1801 * t1808);
t1737 = (-t1797 * t1818 + t1800 * t1810) * t1838 + (t1797 * t1810 + t1800 * t1818) * t1835 + (-t1797 * t1808 + t1800 * t1968) * t1877 + t1874 * (t1797 * t1968 + t1800 * t1808);
t1736 = (-t1796 * t1818 + t1799 * t1809) * t1837 + (t1796 * t1809 + t1799 * t1818) * t1834 + (-t1796 * t1808 + t1799 * t1968) * t1876 + t1873 * (t1796 * t1968 + t1799 * t1808);
t1735 = (-t1842 * t1979 + t1817 + t1904) * t1839 + (t1845 * t1907 + t1804 + t1814) * t1836 + (-t1768 * t1830 - t1767 * t1934 + (t1827 - 0.2e1 * t1958) * t1878 + t1875 * t1824) * pkin(1);
t1734 = (-t1841 * t1980 + t1816 + t1905) * t1838 + (t1844 * t1908 + t1803 + t1813) * t1835 + (-t1766 * t1829 - t1765 * t1935 + (t1826 - 0.2e1 * t1956) * t1877 + t1874 * t1823) * pkin(1);
t1733 = (-t1840 * t1981 + t1815 + t1906) * t1837 + (t1843 * t1909 + t1802 + t1812) * t1834 + (-t1764 * t1828 - t1763 * t1937 + (t1825 - 0.2e1 * t1954) * t1876 + t1873 * t1822) * pkin(1);
t1732 = (t1845 * t1979 + t1804 + t1975) * t1839 + (t1842 * t1907 + t1817 - t1904) * t1836 + (-t1768 * t1934 + t1767 * t1830 + t1824 * t1969 + (-t1957 - t1959) * t1845) * pkin(1);
t1731 = (t1844 * t1980 + t1803 + t1976) * t1838 + (t1841 * t1908 + t1816 - t1905) * t1835 + (-t1766 * t1935 + t1765 * t1829 + t1823 * t1970 + (-t1955 - t1960) * t1844) * pkin(1);
t1730 = (t1843 * t1981 + t1802 + t1977) * t1837 + (t1840 * t1909 + t1815 - t1906) * t1834 + (-t1764 * t1937 + t1763 * t1828 + t1822 * t1971 + (-t1953 - t1961) * t1843) * pkin(1);
t1729 = -t1864 * t1913 - t1941 + (t1753 * t1934 - t1759 * t1830 + t1792 * t1807) * pkin(1) + ((t1899 - t1962) * t1839 + t1898 * t1836) * t1894;
t1728 = -t1864 * t1914 - t1942 + (t1752 * t1935 - t1758 * t1829 + t1788 * t1806) * pkin(1) + ((t1901 - t1963) * t1838 + t1900 * t1835) * t1894;
t1727 = -t1864 * t1915 - t1943 + (t1751 * t1937 - t1757 * t1828 + t1784 * t1805) * pkin(1) + ((t1903 - t1964) * t1837 + t1902 * t1834) * t1894;
t1726 = (t1864 * t1933 - t1931) * qJ(3,1) + (t1753 * t1830 + t1759 * t1934 - t1795 * t1807) * pkin(1) + ((-t1898 + t1827) * t1839 + t1899 * t1836) * t1894;
t1725 = (t1864 * t1930 - t1929) * qJ(3,2) + (t1752 * t1829 + t1758 * t1935 - t1791 * t1806) * pkin(1) + ((-t1900 + t1826) * t1838 + t1901 * t1835) * t1894;
t1724 = (t1864 * t1936 - t1932) * qJ(3,3) + (t1751 * t1828 + t1757 * t1937 - t1787 * t1805) * pkin(1) + ((-t1902 + t1825) * t1837 + t1903 * t1834) * t1894;
t1 = [-m(4) * g(1) + (t1726 * t1747 + t1738 * t1744) * t1750 + (t1725 * t1746 + t1737 * t1743) * t1749 + (t1724 * t1745 + t1736 * t1742) * t1748 + (-t1733 * t1949 - t1734 * t1948 - t1735 * t1947) * m(3); -m(4) * g(2) + (t1729 * t1747 + t1738 * t1741) * t1750 + (t1728 * t1746 + t1737 * t1740) * t1749 + (t1727 * t1745 + t1736 * t1739) * t1748 + (-t1730 * t1949 - t1731 * t1948 - t1732 * t1947) * m(3); m(4) * ((g(1) * t1883 + g(2) * t1882) * t1852 + (g(1) * t1882 - g(2) * t1883) * t1853) + ((t1741 * t1783 + t1744 * t1780) * t1738 + (t1726 * t1780 + t1729 * t1783) * t1747 - (t1732 * t1783 + t1735 * t1780) * m(3) * t1762) * t1750 + ((t1740 * t1782 + t1743 * t1779) * t1737 + (t1725 * t1779 + t1728 * t1782) * t1746 - (t1731 * t1782 + t1734 * t1779) * m(3) * t1761) * t1749 + ((t1739 * t1781 + t1742 * t1778) * t1736 + (t1724 * t1778 + t1727 * t1781) * t1745 - (t1730 * t1781 + t1733 * t1778) * m(3) * t1760) * t1748;];
taugX  = t1;
