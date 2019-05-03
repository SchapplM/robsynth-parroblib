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
%   mass of all robot links (leg links until cut joint, platform)
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
% Datum: 2019-05-03 15:31
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-05-03 15:04:32
% EndTime: 2019-05-03 15:04:34
% DurationCPUTime: 2.60s
% Computational Cost: add. (2328->370), mult. (4797->556), div. (36->3), fcn. (2138->20), ass. (0->252)
t1881 = (qJ(3,3) ^ 2);
t1891 = pkin(1) ^ 2;
t1873 = cos(qJ(1,3));
t1854 = t1873 ^ 2;
t1925 = t1854 * t1891;
t1978 = (2 * t1881) + t1891 - t1925;
t1882 = (qJ(3,2) ^ 2);
t1874 = cos(qJ(1,2));
t1855 = t1874 ^ 2;
t1924 = t1855 * t1891;
t1977 = (2 * t1882) + t1891 - t1924;
t1883 = (qJ(3,1) ^ 2);
t1875 = cos(qJ(1,1));
t1856 = t1875 ^ 2;
t1923 = t1856 * t1891;
t1976 = (2 * t1883) + t1891 - t1923;
t1975 = 2 * pkin(1);
t1867 = legFrame(3,3);
t1837 = sin(t1867);
t1819 = t1837 * qJ(3,3);
t1816 = pkin(2) * t1819;
t1974 = -0.2e1 * t1816;
t1868 = legFrame(2,3);
t1838 = sin(t1868);
t1820 = t1838 * qJ(3,2);
t1817 = pkin(2) * t1820;
t1973 = -0.2e1 * t1817;
t1869 = legFrame(1,3);
t1839 = sin(t1869);
t1821 = t1839 * qJ(3,1);
t1818 = pkin(2) * t1821;
t1972 = -0.2e1 * t1818;
t1870 = sin(qJ(1,3));
t1971 = -0.2e1 * t1870;
t1871 = sin(qJ(1,2));
t1970 = -0.2e1 * t1871;
t1872 = sin(qJ(1,1));
t1969 = -0.2e1 * t1872;
t1968 = -0.2e1 * t1873;
t1967 = -0.2e1 * t1874;
t1966 = -0.2e1 * t1875;
t1965 = m(1) * rSges(1,2);
t1964 = pkin(2) * qJ(3,1);
t1963 = pkin(2) * qJ(3,2);
t1962 = pkin(2) * qJ(3,3);
t1840 = cos(t1867);
t1961 = pkin(2) * t1840;
t1841 = cos(t1868);
t1960 = pkin(2) * t1841;
t1842 = cos(t1869);
t1959 = pkin(2) * t1842;
t1958 = pkin(2) * t1873;
t1957 = pkin(2) * t1874;
t1956 = pkin(2) * t1875;
t1822 = t1837 * pkin(2);
t1823 = t1838 * pkin(2);
t1824 = t1839 * pkin(2);
t1955 = qJ(3,1) * t1842;
t1954 = qJ(3,1) * t1872;
t1953 = qJ(3,2) * t1841;
t1952 = qJ(3,2) * t1871;
t1951 = qJ(3,3) * t1840;
t1950 = qJ(3,3) * t1870;
t1857 = qJ(1,3) + qJ(2,3);
t1831 = sin(t1857);
t1949 = t1831 * qJ(3,3);
t1858 = qJ(1,2) + qJ(2,2);
t1832 = sin(t1858);
t1948 = t1832 * qJ(3,2);
t1859 = qJ(1,1) + qJ(2,1);
t1833 = sin(t1859);
t1947 = t1833 * qJ(3,1);
t1834 = cos(t1857);
t1825 = t1834 ^ 2;
t1890 = pkin(2) ^ 2;
t1843 = -t1881 + t1890;
t1828 = 0.1e1 + t1843;
t1860 = 0.2e1 + t1891;
t1922 = t1870 * t1873;
t1892 = -t1828 * t1854 + 0.2e1 * t1922 * t1962;
t1919 = -0.1e1 / 0.2e1 - t1890 / 0.2e1;
t1745 = 0.1e1 / ((pkin(1) * (t1828 * t1922 + (0.2e1 * t1854 - 0.1e1) * t1962) * t1831 + t1950) * t1834 * t1975 + pkin(1) * t1949 * t1968 - t1881 * t1860 + (0.2e1 * (t1881 / 0.2e1 - t1892 + t1919) * t1825 + t1892) * t1891);
t1793 = -g(1) * t1837 + g(2) * t1840;
t1796 = g(1) * t1840 + g(2) * t1837;
t1757 = -t1793 * t1834 + t1796 * t1831;
t1946 = t1745 * t1757;
t1835 = cos(t1858);
t1826 = t1835 ^ 2;
t1844 = -t1882 + t1890;
t1829 = 0.1e1 + t1844;
t1921 = t1871 * t1874;
t1893 = -t1829 * t1855 + 0.2e1 * t1921 * t1963;
t1746 = 0.1e1 / ((pkin(1) * (t1829 * t1921 + (0.2e1 * t1855 - 0.1e1) * t1963) * t1832 + t1952) * t1835 * t1975 + pkin(1) * t1948 * t1967 - t1882 * t1860 + (0.2e1 * (t1882 / 0.2e1 - t1893 + t1919) * t1826 + t1893) * t1891);
t1794 = -g(1) * t1838 + g(2) * t1841;
t1797 = g(1) * t1841 + g(2) * t1838;
t1758 = -t1794 * t1835 + t1797 * t1832;
t1945 = t1746 * t1758;
t1836 = cos(t1859);
t1827 = t1836 ^ 2;
t1845 = -t1883 + t1890;
t1830 = 0.1e1 + t1845;
t1920 = t1872 * t1875;
t1894 = -t1830 * t1856 + 0.2e1 * t1920 * t1964;
t1747 = 0.1e1 / ((pkin(1) * (t1830 * t1920 + (0.2e1 * t1856 - 0.1e1) * t1964) * t1833 + t1954) * t1836 * t1975 + pkin(1) * t1947 * t1966 - t1883 * t1860 + (0.2e1 * (t1883 / 0.2e1 - t1894 + t1919) * t1827 + t1894) * t1891);
t1795 = -g(1) * t1839 + g(2) * t1842;
t1798 = g(1) * t1842 + g(2) * t1839;
t1759 = -t1795 * t1836 + t1798 * t1833;
t1944 = t1747 * t1759;
t1781 = t1819 + t1961;
t1943 = t1781 * t1870;
t1785 = t1820 + t1960;
t1942 = t1785 * t1871;
t1789 = t1821 + t1959;
t1941 = t1789 * t1872;
t1940 = t1819 * t1834;
t1939 = t1820 * t1835;
t1938 = t1821 * t1836;
t1937 = t1828 * t1840;
t1936 = t1829 * t1841;
t1935 = t1830 * t1842;
t1934 = t1831 * t1834;
t1933 = t1832 * t1835;
t1932 = t1833 * t1836;
t1931 = t1833 * t1839;
t1930 = t1836 * t1842;
t1929 = t1837 * t1831;
t1928 = t1838 * t1832;
t1927 = t1840 * t1834;
t1926 = t1841 * t1835;
t1865 = t1890 / 0.2e1;
t1918 = t1865 - t1881 / 0.2e1;
t1917 = t1865 - t1882 / 0.2e1;
t1916 = t1865 - t1883 / 0.2e1;
t1915 = pkin(2) * t1955;
t1914 = pkin(2) * t1953;
t1913 = pkin(2) * t1951;
t1912 = t1840 * t1949;
t1911 = t1841 * t1948;
t1910 = t1842 * t1947;
t1909 = t1891 * t1922;
t1908 = t1891 * t1921;
t1907 = t1891 * t1920;
t1906 = -(2 * t1881) - t1925;
t1905 = -(2 * t1882) - t1924;
t1904 = -(2 * t1883) - t1923;
t1903 = t1840 * t1909;
t1902 = t1841 * t1908;
t1901 = t1842 * t1907;
t1782 = t1819 - t1961;
t1783 = t1822 + t1951;
t1900 = -t1782 * t1854 - t1783 * t1922;
t1899 = -t1782 * t1922 + t1783 * t1854;
t1786 = t1820 - t1960;
t1787 = t1823 + t1953;
t1898 = -t1786 * t1855 - t1787 * t1921;
t1897 = -t1786 * t1921 + t1787 * t1855;
t1790 = t1821 - t1959;
t1791 = t1824 + t1955;
t1896 = -t1790 * t1856 - t1791 * t1920;
t1895 = -t1790 * t1920 + t1791 * t1856;
t1889 = koppelP(1,1);
t1888 = koppelP(2,1);
t1887 = koppelP(3,1);
t1886 = koppelP(1,2);
t1885 = koppelP(2,2);
t1884 = koppelP(3,2);
t1880 = rSges(4,1);
t1879 = rSges(4,2);
t1878 = xP(3);
t1877 = m(2) * rSges(2,2);
t1876 = pkin(2) + rSges(3,1);
t1864 = qJ(3,1) + rSges(3,3);
t1863 = qJ(3,2) + rSges(3,3);
t1862 = qJ(3,3) + rSges(3,3);
t1861 = 0.1e1 + t1891;
t1850 = cos(t1878);
t1849 = sin(t1878);
t1815 = m(2) * rSges(2,1) + m(3) * t1876;
t1814 = -0.2e1 * t1915;
t1813 = -0.2e1 * t1914;
t1812 = -0.2e1 * t1913;
t1811 = 0.2e1 * t1818;
t1810 = 0.2e1 * t1817;
t1809 = 0.2e1 * t1816;
t1808 = -m(3) * t1864 + t1877;
t1807 = -m(3) * t1863 + t1877;
t1806 = -m(3) * t1862 + t1877;
t1805 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1804 = -t1954 + t1956;
t1803 = -t1952 + t1957;
t1802 = -t1950 + t1958;
t1801 = t1839 * t1907;
t1800 = t1838 * t1908;
t1799 = t1837 * t1909;
t1792 = t1824 - t1955;
t1788 = t1823 - t1953;
t1784 = t1822 - t1951;
t1780 = -t1849 * t1886 + t1850 * t1889;
t1779 = -t1849 * t1885 + t1850 * t1888;
t1778 = -t1849 * t1884 + t1850 * t1887;
t1777 = -t1849 * t1889 - t1850 * t1886;
t1776 = -t1849 * t1888 - t1850 * t1885;
t1775 = -t1849 * t1887 - t1850 * t1884;
t1774 = t1839 * t1845 + t1814;
t1773 = t1838 * t1844 + t1813;
t1772 = t1837 * t1843 + t1812;
t1771 = t1830 * t1839 + t1814;
t1770 = t1829 * t1838 + t1813;
t1769 = t1828 * t1837 + t1812;
t1768 = t1771 * t1875;
t1767 = t1770 * t1874;
t1766 = t1769 * t1873;
t1765 = t1792 * t1875 - t1941;
t1764 = t1789 * t1875 + t1792 * t1872;
t1763 = t1788 * t1874 - t1942;
t1762 = t1785 * t1874 + t1788 * t1871;
t1761 = t1784 * t1873 - t1943;
t1760 = t1781 * t1873 + t1784 * t1870;
t1756 = (t1842 * t1845 + t1811) * t1875 + t1774 * t1872;
t1755 = (t1841 * t1844 + t1810) * t1874 + t1773 * t1871;
t1754 = (t1840 * t1843 + t1809) * t1873 + t1772 * t1870;
t1753 = (t1811 + t1935) * t1875 + t1872 * t1771;
t1752 = (t1810 + t1936) * t1874 + t1871 * t1770;
t1751 = (t1809 + t1937) * t1873 + t1870 * t1769;
t1750 = t1774 * t1875 + (t1842 * t1916 + t1818) * t1969;
t1749 = t1773 * t1874 + (t1841 * t1917 + t1817) * t1970;
t1748 = t1772 * t1873 + (t1840 * t1918 + t1816) * t1971;
t1744 = ((-t1795 * t1876 - t1798 * t1864) * m(3) + m(2) * (-rSges(2,1) * t1795 + rSges(2,2) * t1798)) * t1836 + t1833 * ((-t1795 * t1864 + t1798 * t1876) * m(3) + m(2) * (rSges(2,1) * t1798 + rSges(2,2) * t1795));
t1743 = ((-t1794 * t1876 - t1797 * t1863) * m(3) + m(2) * (-rSges(2,1) * t1794 + rSges(2,2) * t1797)) * t1835 + t1832 * ((-t1794 * t1863 + t1797 * t1876) * m(3) + m(2) * (rSges(2,1) * t1797 + rSges(2,2) * t1794));
t1742 = ((-t1793 * t1876 - t1796 * t1862) * m(3) + m(2) * (-rSges(2,1) * t1793 + rSges(2,2) * t1796)) * t1834 + t1831 * ((-t1793 * t1862 + t1796 * t1876) * m(3) + m(2) * (rSges(2,1) * t1796 + rSges(2,2) * t1793));
t1741 = (-t1930 + t1931) * qJ(3,1) + (-t1768 * t1827 - t1753 * t1932 + (t1839 * t1890 + t1839 - t1915) * t1875 + (-(t1972 - t1935) * t1827 - qJ(3,1) * t1792) * t1872) * pkin(1);
t1740 = (-t1926 + t1928) * qJ(3,2) + (-t1767 * t1826 - t1752 * t1933 + (t1838 * t1890 + t1838 - t1914) * t1874 + (-(t1973 - t1936) * t1826 - qJ(3,2) * t1788) * t1871) * pkin(1);
t1739 = (-t1927 + t1929) * qJ(3,3) + (-t1766 * t1825 - t1751 * t1934 + (t1837 * t1890 + t1837 - t1913) * t1873 + (-(t1974 - t1937) * t1825 - qJ(3,3) * t1784) * t1870) * pkin(1);
t1738 = -t1910 - t1938 + (t1753 * t1827 - (t1768 + ((0.1e1 / 0.2e1 + t1916) * t1842 + t1818) * t1969) * t1932 - (t1842 * t1890 + t1818 + t1842) * t1875 + qJ(3,1) * t1941) * pkin(1);
t1737 = -t1911 - t1939 + (t1752 * t1826 - (t1767 + ((0.1e1 / 0.2e1 + t1917) * t1841 + t1817) * t1970) * t1933 - (t1841 * t1890 + t1817 + t1841) * t1874 + qJ(3,2) * t1942) * pkin(1);
t1736 = -t1912 - t1940 + (t1751 * t1825 - (t1766 + ((0.1e1 / 0.2e1 + t1918) * t1840 + t1816) * t1971) * t1934 - (t1840 * t1890 + t1816 + t1840) * t1873 + qJ(3,3) * t1943) * pkin(1);
t1735 = (-t1795 * t1815 + t1798 * t1808) * t1836 + (t1795 * t1808 + t1798 * t1815) * t1833 + (-t1795 * t1805 + t1798 * t1965) * t1875 + (t1795 * t1965 + t1798 * t1805) * t1872;
t1734 = (-t1794 * t1815 + t1797 * t1807) * t1835 + (t1794 * t1807 + t1797 * t1815) * t1832 + (-t1794 * t1805 + t1797 * t1965) * t1874 + (t1794 * t1965 + t1797 * t1805) * t1871;
t1733 = (-t1793 * t1815 + t1796 * t1806) * t1834 + (t1793 * t1806 + t1796 * t1815) * t1831 + (-t1793 * t1805 + t1796 * t1965) * t1873 + (t1793 * t1965 + t1796 * t1805) * t1870;
t1732 = (-t1839 * t1976 + t1814 + t1901) * t1836 + (t1842 * t1904 + t1801 + t1811) * t1833 + (-t1765 * t1827 - t1764 * t1932 + (t1824 - 0.2e1 * t1955) * t1875 + t1872 * t1821) * pkin(1);
t1731 = (-t1838 * t1977 + t1813 + t1902) * t1835 + (t1841 * t1905 + t1800 + t1810) * t1832 + (-t1763 * t1826 - t1762 * t1933 + (t1823 - 0.2e1 * t1953) * t1874 + t1871 * t1820) * pkin(1);
t1730 = (-t1837 * t1978 + t1812 + t1903) * t1834 + (t1840 * t1906 + t1799 + t1809) * t1831 + (-t1761 * t1825 - t1760 * t1934 + (t1822 - 0.2e1 * t1951) * t1873 + t1870 * t1819) * pkin(1);
t1729 = (t1842 * t1976 + t1801 + t1972) * t1836 + (t1839 * t1904 + t1814 - t1901) * t1833 + (-t1765 * t1932 + t1764 * t1827 + t1821 * t1966 + (-t1954 - t1956) * t1842) * pkin(1);
t1728 = (t1841 * t1977 + t1800 + t1973) * t1835 + (t1838 * t1905 + t1813 - t1902) * t1832 + (-t1763 * t1933 + t1762 * t1826 + t1820 * t1967 + (-t1952 - t1957) * t1841) * pkin(1);
t1727 = (t1840 * t1978 + t1799 + t1974) * t1834 + (t1837 * t1906 + t1812 - t1903) * t1831 + (-t1761 * t1934 + t1760 * t1825 + t1819 * t1968 + (-t1950 - t1958) * t1840) * pkin(1);
t1726 = -t1861 * t1910 - t1938 + (t1750 * t1932 - t1756 * t1827 + t1789 * t1804) * pkin(1) + ((t1896 - t1959) * t1836 + t1895 * t1833) * t1891;
t1725 = -t1861 * t1911 - t1939 + (t1749 * t1933 - t1755 * t1826 + t1785 * t1803) * pkin(1) + ((t1898 - t1960) * t1835 + t1897 * t1832) * t1891;
t1724 = -t1861 * t1912 - t1940 + (t1748 * t1934 - t1754 * t1825 + t1781 * t1802) * pkin(1) + ((t1900 - t1961) * t1834 + t1899 * t1831) * t1891;
t1723 = (t1861 * t1931 - t1930) * qJ(3,1) + (t1750 * t1827 + t1756 * t1932 - t1792 * t1804) * pkin(1) + ((-t1895 + t1824) * t1836 + t1896 * t1833) * t1891;
t1722 = (t1861 * t1928 - t1926) * qJ(3,2) + (t1749 * t1826 + t1755 * t1933 - t1788 * t1803) * pkin(1) + ((-t1897 + t1823) * t1835 + t1898 * t1832) * t1891;
t1721 = (t1861 * t1929 - t1927) * qJ(3,3) + (t1748 * t1825 + t1754 * t1934 - t1784 * t1802) * pkin(1) + ((-t1899 + t1822) * t1834 + t1900 * t1831) * t1891;
t1 = [-m(4) * g(1) + (t1723 * t1744 + t1735 * t1741) * t1747 + (t1722 * t1743 + t1734 * t1740) * t1746 + (t1721 * t1742 + t1733 * t1739) * t1745 + (-t1730 * t1946 - t1731 * t1945 - t1732 * t1944) * m(3); -m(4) * g(2) + (t1726 * t1744 + t1735 * t1738) * t1747 + (t1725 * t1743 + t1734 * t1737) * t1746 + (t1724 * t1742 + t1733 * t1736) * t1745 + (-t1727 * t1946 - t1728 * t1945 - t1729 * t1944) * m(3); ((g(1) * t1880 + g(2) * t1879) * t1849 + (g(1) * t1879 - g(2) * t1880) * t1850) * m(4) + ((t1738 * t1780 + t1741 * t1777) * t1735 + (t1723 * t1777 + t1726 * t1780) * t1744 - (t1729 * t1780 + t1732 * t1777) * m(3) * t1759) * t1747 + ((t1737 * t1779 + t1740 * t1776) * t1734 + (t1722 * t1776 + t1725 * t1779) * t1743 - (t1728 * t1779 + t1731 * t1776) * m(3) * t1758) * t1746 + ((t1736 * t1778 + t1739 * t1775) * t1733 + (t1721 * t1775 + t1724 * t1778) * t1742 - (t1727 * t1778 + t1730 * t1775) * m(3) * t1757) * t1745;];
taugX  = t1;
