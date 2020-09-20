% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR9V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR9V1G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:47:51
% EndTime: 2020-08-06 18:48:07
% DurationCPUTime: 16.32s
% Computational Cost: add. (111009->747), mult. (110965->1122), div. (9828->14), fcn. (53133->107), ass. (0->522)
t1910 = (pkin(7) + qJ(3,1));
t1848 = 2 * t1910;
t1835 = sin(t1848);
t1946 = 2 * pkin(7);
t1907 = t1946 + qJ(3,1);
t1866 = sin(t1907);
t1930 = sin(qJ(3,1));
t2109 = t1866 + t1930;
t1869 = sin(t1910);
t2241 = 0.2e1 * t1869;
t1751 = pkin(1) * t2241 + pkin(2) * t2109 + pkin(3) * t1835;
t1875 = cos(t1910);
t1861 = 0.1e1 / t1875;
t1923 = xDDP(3);
t2064 = t1751 * t1861 * t1923;
t1916 = cos(pkin(7));
t1845 = t1916 * pkin(2) + pkin(1);
t1922 = pkin(5) + qJ(2,1);
t1900 = -pkin(6) - t1922;
t1931 = sin(qJ(1,1));
t1937 = cos(qJ(1,1));
t1763 = t1845 * t1931 + t1900 * t1937;
t1766 = t1845 * t1937 - t1900 * t1931;
t1919 = legFrame(1,3);
t1883 = sin(t1919);
t1886 = cos(t1919);
t1775 = t1883 * t1937 + t1886 * t1931;
t2203 = pkin(3) * t1875;
t1727 = t1763 * t1886 + t1766 * t1883 + t1775 * t2203;
t1924 = xDDP(2);
t2159 = t1727 * t1924;
t2006 = t1883 * t1931 - t1886 * t1937;
t1724 = -t1763 * t1883 + t1766 * t1886 - t2006 * t2203;
t1925 = xDDP(1);
t2162 = t1724 * t1925;
t2272 = t2064 / 0.2e1 + t2159 + t2162;
t1909 = pkin(7) + qJ(3,2);
t1847 = 2 * t1909;
t1834 = sin(t1847);
t1906 = t1946 + qJ(3,2);
t1865 = sin(t1906);
t1928 = sin(qJ(3,2));
t2110 = t1865 + t1928;
t1868 = sin(t1909);
t2242 = 0.2e1 * t1868;
t1750 = pkin(1) * t2242 + pkin(2) * t2110 + pkin(3) * t1834;
t1874 = cos(t1909);
t1858 = 0.1e1 / t1874;
t2065 = t1750 * t1858 * t1923;
t1921 = pkin(5) + qJ(2,2);
t1899 = -pkin(6) - t1921;
t1929 = sin(qJ(1,2));
t1935 = cos(qJ(1,2));
t1762 = t1845 * t1929 + t1899 * t1935;
t1765 = t1845 * t1935 - t1899 * t1929;
t1918 = legFrame(2,3);
t1882 = sin(t1918);
t1885 = cos(t1918);
t1774 = t1882 * t1935 + t1885 * t1929;
t2204 = pkin(3) * t1874;
t1726 = t1762 * t1885 + t1765 * t1882 + t1774 * t2204;
t2160 = t1726 * t1924;
t2007 = t1882 * t1929 - t1885 * t1935;
t1723 = -t1762 * t1882 + t1765 * t1885 - t2007 * t2204;
t2163 = t1723 * t1925;
t2271 = t2065 / 0.2e1 + t2160 + t2163;
t1908 = pkin(7) + qJ(3,3);
t1846 = 2 * t1908;
t1833 = sin(t1846);
t1905 = t1946 + qJ(3,3);
t1864 = sin(t1905);
t1926 = sin(qJ(3,3));
t2111 = t1864 + t1926;
t1867 = sin(t1908);
t2243 = 0.2e1 * t1867;
t1749 = pkin(1) * t2243 + pkin(2) * t2111 + pkin(3) * t1833;
t1873 = cos(t1908);
t1855 = 0.1e1 / t1873;
t2066 = t1749 * t1855 * t1923;
t1920 = pkin(5) + qJ(2,3);
t1898 = -pkin(6) - t1920;
t1927 = sin(qJ(1,3));
t1933 = cos(qJ(1,3));
t1761 = t1845 * t1927 + t1898 * t1933;
t1764 = t1845 * t1933 - t1898 * t1927;
t1917 = legFrame(3,3);
t1881 = sin(t1917);
t1884 = cos(t1917);
t1773 = t1881 * t1933 + t1884 * t1927;
t2205 = pkin(3) * t1873;
t1725 = t1761 * t1884 + t1764 * t1881 + t1773 * t2205;
t2161 = t1725 * t1924;
t2008 = t1881 * t1927 - t1884 * t1933;
t1722 = -t1761 * t1881 + t1764 * t1884 - t2008 * t2205;
t2164 = t1722 * t1925;
t2270 = t2066 / 0.2e1 + t2161 + t2164;
t1836 = cos(t1846);
t1932 = cos(qJ(3,3));
t1938 = xDP(3);
t1939 = xDP(2);
t1940 = xDP(1);
t1731 = ((t1927 * t1939 + t1933 * t1940) * t1884 + t1881 * (-t1927 * t1940 + t1933 * t1939)) * t1873 + t1867 * t1938;
t1888 = 0.1e1 / t1898;
t2137 = t1855 * t1888;
t2072 = t1731 * t2137;
t2038 = pkin(3) * t2072;
t2062 = t1873 * t2137;
t2202 = t1731 * pkin(1);
t1870 = cos(t1905);
t2208 = pkin(2) * t1870;
t1876 = t1917 + qJ(1,3);
t1841 = qJ(3,3) + t1876;
t1821 = pkin(7) + t1841;
t1803 = sin(t1821);
t1842 = -qJ(3,3) + t1876;
t1822 = -pkin(7) + t1842;
t1804 = sin(t1822);
t1809 = cos(t1821);
t1810 = cos(t1822);
t1948 = 2 * qJ(3,3);
t2089 = t1946 + t1876;
t1815 = t1948 + t2089;
t1816 = qJ(3,3) + t2089;
t1947 = -2 * pkin(7);
t2088 = t1947 + t1876;
t1817 = -qJ(3,3) + t2088;
t1818 = -(2 * qJ(3,3)) + t2088;
t2250 = 0.2e1 * pkin(1);
t1879 = t1939 * t2250;
t1897 = pkin(1) * t1940;
t2266 = 0.4e1 * pkin(1);
t2104 = t1938 * t2266;
t2235 = -0.2e1 * t1940;
t2263 = 0.2e1 * t1938;
t1707 = t1867 * t2104 + (t1803 + t1804) * (t1898 * t2235 + t1879) + 0.2e1 * (t1809 + t1810) * (t1898 * t1939 + t1897) + (t1833 * t2263 + (cos(t1818) + cos(t1815) + 0.2e1 * cos(t1876)) * t1940 + (sin(t1818) + sin(t1815) + 0.2e1 * sin(t1876)) * t1939) * pkin(3) + (t2111 * t2263 + (cos(t1817) + cos(t1816) + cos(t1842) + cos(t1841)) * t1940 + (sin(t1817) + sin(t1816) + sin(t1842) + sin(t1841)) * t1939) * pkin(2);
t2228 = t1707 / 0.2e1;
t1695 = t1836 * t2038 - (-0.2e1 * t2202 + t2228) * t2062 + (t1932 * pkin(2) + pkin(3) + t2208) * t2072;
t1755 = t2008 * t1888 * t1925;
t1758 = t1773 * t1888 * t1924;
t2063 = t1867 * t2137;
t1767 = t1923 * t2063;
t1856 = 0.1e1 / t1873 ^ 2;
t1857 = t1855 * t1856;
t1914 = t1938 ^ 2;
t1952 = 0.1e1 / pkin(3);
t2125 = t1914 * t1952;
t1779 = t1857 * t1888 * t2125;
t1961 = t1898 ^ 2;
t1889 = 0.1e1 / t1961;
t1890 = t1888 * t1889;
t2155 = t1731 * t1856;
t2049 = t2155 / 0.4e1;
t2135 = t1856 * t1889;
t2071 = t1731 * t2135;
t1686 = t1755 - t1758 - t1767 - t1695 * t2071 / 0.2e1 + t1890 * t1707 * t2049 - t1779;
t1683 = t1686 * pkin(1);
t2118 = t1952 * t1938;
t2167 = (-t1731 * t1833 / 0.2e1 + (t1845 + t2205) * t2118) * t1855 * t1938;
t2083 = t1888 * t2167;
t2124 = t1914 / pkin(3) ^ 2;
t2096 = -0.2e1 * t2124;
t2259 = (t1870 + t1932) * pkin(2);
t1705 = (t2250 + ((t1836 + 0.1e1) * pkin(3) + t2259) * t1855) * t1888 * t1731;
t2169 = t1705 * t1707;
t1945 = 3 * pkin(7);
t1951 = pkin(3) ^ 2;
t1954 = pkin(2) ^ 2;
t2023 = pkin(2) * t2038;
t2042 = 0.2e1 * pkin(3) * t2118;
t2043 = -0.4e1 * pkin(1) ^ 2 - 0.3e1 * t1951 - 0.2e1 * t1954;
t2246 = 0.4e1 * (t2202 - t1707 / 0.8e1) * t2137;
t2100 = pkin(3) * t2246;
t2240 = 0.4e1 * t1916;
t1668 = -t1898 * t1855 * t1833 * t2042 + t1836 * t2100 - ((-0.4e1 * t1961 + t2043) * t1731 + pkin(1) * t1707) * t2062 + t2023 * t2240 + t2100 + t2246 * t2259 + 0.2e1 * (cos((t1948 + t1945)) + cos((t1948 + pkin(7)))) * t2023 + (t1951 * cos((3 * t1908)) + (cos((qJ(3,3) - pkin(7))) + cos((t1945 + qJ(3,3)))) * t1954) * t2072;
t2182 = t1668 * t1731;
t2269 = (t1889 * (t2182 / 0.2e1 - t2169 / 0.4e1) + t1920 * t2096 + 0.2e1 * t2083) * t1856 + (t2066 + 0.2e1 * t2161 + 0.2e1 * t2164) * t1888 + 0.4e1 * t1683;
t1837 = cos(t1847);
t1934 = cos(qJ(3,2));
t1732 = ((t1929 * t1939 + t1935 * t1940) * t1885 + t1882 * (-t1929 * t1940 + t1935 * t1939)) * t1874 + t1868 * t1938;
t1891 = 0.1e1 / t1899;
t2134 = t1858 * t1891;
t2070 = t1732 * t2134;
t2037 = pkin(3) * t2070;
t2060 = t1874 * t2134;
t2201 = t1732 * pkin(1);
t1871 = cos(t1906);
t2207 = pkin(2) * t1871;
t1877 = t1918 + qJ(1,2);
t1843 = qJ(3,2) + t1877;
t1825 = pkin(7) + t1843;
t1805 = sin(t1825);
t1844 = -qJ(3,2) + t1877;
t1826 = -pkin(7) + t1844;
t1806 = sin(t1826);
t1811 = cos(t1825);
t1812 = cos(t1826);
t1949 = 2 * qJ(3,2);
t2091 = t1946 + t1877;
t1823 = t1949 + t2091;
t2090 = t1947 + t1877;
t1824 = -(2 * qJ(3,2)) + t2090;
t1827 = qJ(3,2) + t2091;
t1828 = -qJ(3,2) + t2090;
t1708 = t1868 * t2104 + (t1805 + t1806) * (t1899 * t2235 + t1879) + 0.2e1 * (t1811 + t1812) * (t1899 * t1939 + t1897) + (t1834 * t2263 + (cos(t1824) + cos(t1823) + 0.2e1 * cos(t1877)) * t1940 + (sin(t1824) + sin(t1823) + 0.2e1 * sin(t1877)) * t1939) * pkin(3) + (t2110 * t2263 + (cos(t1828) + cos(t1827) + cos(t1844) + cos(t1843)) * t1940 + (sin(t1828) + sin(t1827) + sin(t1844) + sin(t1843)) * t1939) * pkin(2);
t2227 = t1708 / 0.2e1;
t1696 = t1837 * t2037 - (-0.2e1 * t2201 + t2227) * t2060 + (pkin(2) * t1934 + pkin(3) + t2207) * t2070;
t1756 = t2007 * t1891 * t1925;
t1759 = t1774 * t1891 * t1924;
t2061 = t1868 * t2134;
t1768 = t1923 * t2061;
t1859 = 0.1e1 / t1874 ^ 2;
t1860 = t1858 * t1859;
t1780 = t1860 * t1891 * t2125;
t1963 = t1899 ^ 2;
t1892 = 0.1e1 / t1963;
t1893 = t1891 * t1892;
t2154 = t1732 * t1859;
t2048 = t2154 / 0.4e1;
t2132 = t1859 * t1892;
t2069 = t1732 * t2132;
t1687 = t1756 - t1759 - t1768 - t1696 * t2069 / 0.2e1 + t1893 * t1708 * t2048 - t1780;
t1684 = t1687 * pkin(1);
t2166 = (-t1732 * t1834 / 0.2e1 + (t1845 + t2204) * t2118) * t1858 * t1938;
t2081 = t1891 * t2166;
t2258 = (t1871 + t1934) * pkin(2);
t1706 = (t2250 + ((t1837 + 0.1e1) * pkin(3) + t2258) * t1858) * t1891 * t1732;
t2168 = t1706 * t1708;
t2022 = pkin(2) * t2037;
t2245 = 0.4e1 * (t2201 - t1708 / 0.8e1) * t2134;
t2099 = pkin(3) * t2245;
t1669 = -t1899 * t1858 * t1834 * t2042 + t1837 * t2099 - ((-0.4e1 * t1963 + t2043) * t1732 + pkin(1) * t1708) * t2060 + t2022 * t2240 + t2099 + t2245 * t2258 + 0.2e1 * (cos((t1949 + t1945)) + cos((t1949 + pkin(7)))) * t2022 + (t1951 * cos((3 * t1909)) + (cos((qJ(3,2) - pkin(7))) + cos((t1945 + qJ(3,2)))) * t1954) * t2070;
t2181 = t1669 * t1732;
t2268 = (t1892 * (t2181 / 0.2e1 - t2168 / 0.4e1) + t1921 * t2096 + 0.2e1 * t2081) * t1859 + (t2065 + 0.2e1 * t2160 + 0.2e1 * t2163) * t1891 + 0.4e1 * t1684;
t1838 = cos(t1848);
t1936 = cos(qJ(3,1));
t1733 = ((t1931 * t1939 + t1937 * t1940) * t1886 + t1883 * (-t1931 * t1940 + t1937 * t1939)) * t1875 + t1869 * t1938;
t1894 = 0.1e1 / t1900;
t2131 = t1861 * t1894;
t2068 = t1733 * t2131;
t2036 = pkin(3) * t2068;
t2058 = t1875 * t2131;
t2200 = t1733 * pkin(1);
t1872 = cos(t1907);
t2206 = pkin(2) * t1872;
t1878 = t1919 + qJ(1,1);
t1839 = qJ(3,1) + t1878;
t1829 = pkin(7) + t1839;
t1807 = sin(t1829);
t1840 = -qJ(3,1) + t1878;
t1830 = -pkin(7) + t1840;
t1808 = sin(t1830);
t1813 = cos(t1829);
t1814 = cos(t1830);
t2093 = t1946 + t1878;
t1819 = qJ(3,1) + t2093;
t2092 = t1947 + t1878;
t1820 = -qJ(3,1) + t2092;
t1950 = 2 * qJ(3,1);
t1831 = t1950 + t2093;
t1832 = -(2 * qJ(3,1)) + t2092;
t1709 = t1869 * t2104 + (t1807 + t1808) * (t1900 * t2235 + t1879) + 0.2e1 * (t1813 + t1814) * (t1900 * t1939 + t1897) + (t1835 * t2263 + (cos(t1831) + 0.2e1 * cos(t1878) + cos(t1832)) * t1940 + (sin(t1832) + sin(t1831) + 0.2e1 * sin(t1878)) * t1939) * pkin(3) + (t2109 * t2263 + (cos(t1820) + cos(t1819) + cos(t1840) + cos(t1839)) * t1940 + (sin(t1820) + sin(t1819) + sin(t1840) + sin(t1839)) * t1939) * pkin(2);
t2226 = t1709 / 0.2e1;
t1697 = t1838 * t2036 - (-0.2e1 * t2200 + t2226) * t2058 + (pkin(2) * t1936 + pkin(3) + t2206) * t2068;
t1757 = t2006 * t1894 * t1925;
t1760 = t1775 * t1894 * t1924;
t2059 = t1869 * t2131;
t1769 = t1923 * t2059;
t1862 = 0.1e1 / t1875 ^ 2;
t1863 = t1861 * t1862;
t1781 = t1863 * t1894 * t2125;
t1965 = t1900 ^ 2;
t1895 = 0.1e1 / t1965;
t1896 = t1894 * t1895;
t2153 = t1733 * t1862;
t2047 = t2153 / 0.4e1;
t2129 = t1862 * t1895;
t2067 = t1733 * t2129;
t1688 = t1757 - t1760 - t1769 - t1697 * t2067 / 0.2e1 + t1896 * t1709 * t2047 - t1781;
t1685 = t1688 * pkin(1);
t2165 = (-t1733 * t1835 / 0.2e1 + (t1845 + t2203) * t2118) * t1861 * t1938;
t2079 = t1894 * t2165;
t2257 = (t1872 + t1936) * pkin(2);
t1704 = (t2250 + ((t1838 + 0.1e1) * pkin(3) + t2257) * t1861) * t1894 * t1733;
t2170 = t1704 * t1709;
t2021 = pkin(2) * t2036;
t2244 = 0.4e1 * (t2200 - t1709 / 0.8e1) * t2131;
t2098 = pkin(3) * t2244;
t1670 = -t1900 * t1861 * t1835 * t2042 + t1838 * t2098 - ((-0.4e1 * t1965 + t2043) * t1733 + pkin(1) * t1709) * t2058 + t2021 * t2240 + t2098 + t2244 * t2257 + 0.2e1 * (cos((t1945 + t1950)) + cos((pkin(7) + t1950))) * t2021 + (t1951 * cos((3 * t1910)) + (cos((qJ(3,1) - pkin(7))) + cos((qJ(3,1) + t1945))) * t1954) * t2068;
t2180 = t1670 * t1733;
t2267 = (t1895 * (t2180 / 0.2e1 - t2170 / 0.4e1) + t1922 * t2096 + 0.2e1 * t2079) * t1862 + (t2064 + 0.2e1 * t2159 + 0.2e1 * t2162) * t1894 + 0.4e1 * t1685;
t1915 = sin(pkin(7));
t2265 = 0.2e1 * t1915;
t2264 = -0.2e1 * t1916;
t2080 = t1862 * t2165;
t2262 = t1670 * t1895 * t2047 - t2129 * t2170 / 0.8e1 + (t2080 + t2272) * t1894;
t2082 = t1859 * t2166;
t2261 = t1669 * t1892 * t2048 - t2132 * t2168 / 0.8e1 + (t2082 + t2271) * t1891;
t2084 = t1856 * t2167;
t2260 = t1668 * t1889 * t2049 - t2135 * t2169 / 0.8e1 + (t2084 + t2270) * t1888;
t2188 = t1886 * g(2);
t2195 = t1883 * g(1);
t1784 = -t2188 + t2195;
t2189 = t1886 * g(1);
t2194 = t1883 * g(2);
t1787 = t2189 + t2194;
t1748 = -t1784 * t1931 + t1787 * t1937;
t1745 = t1784 * t1937 + t1787 * t1931;
t2190 = t1885 * g(2);
t2197 = t1882 * g(1);
t1783 = -t2190 + t2197;
t2191 = t1885 * g(1);
t2196 = t1882 * g(2);
t1786 = t2191 + t2196;
t1747 = -t1783 * t1929 + t1786 * t1935;
t1744 = t1783 * t1935 + t1786 * t1929;
t2192 = t1884 * g(2);
t2199 = t1881 * g(1);
t1782 = -t2192 + t2199;
t2193 = t1884 * g(1);
t2198 = t1881 * g(2);
t1785 = t2193 + t2198;
t1746 = -t1782 * t1927 + t1785 * t1933;
t1743 = t1782 * t1933 + t1785 * t1927;
t2087 = t1707 / 0.16e2;
t1659 = (t2193 / 0.2e1 + t2198 / 0.2e1) * t1927 + t1683 + (-t2192 / 0.2e1 + t2199 / 0.2e1) * t1933 + (t2182 / 0.8e1 - t1705 * t2087) * t2135 - (-t2164 / 0.2e1 - t2161 / 0.2e1 - t2066 / 0.4e1 - t2084 / 0.2e1) * t1888;
t2249 = 0.2e1 * t1659;
t2086 = t1708 / 0.16e2;
t1660 = (t2191 / 0.2e1 + t2196 / 0.2e1) * t1929 + t1684 + (-t2190 / 0.2e1 + t2197 / 0.2e1) * t1935 + (t2181 / 0.8e1 - t1706 * t2086) * t2132 - (-t2163 / 0.2e1 - t2160 / 0.2e1 - t2065 / 0.4e1 - t2082 / 0.2e1) * t1891;
t2248 = 0.2e1 * t1660;
t2085 = t1709 / 0.16e2;
t1661 = (t2189 / 0.2e1 + t2194 / 0.2e1) * t1931 + t1685 + (-t2188 / 0.2e1 + t2195 / 0.2e1) * t1937 + (t2180 / 0.8e1 - t1704 * t2085) * t2129 - (-t2162 / 0.2e1 - t2159 / 0.2e1 - t2064 / 0.4e1 - t2080 / 0.2e1) * t1894;
t2247 = 0.2e1 * t1661;
t2239 = -0.2e1 * t1920;
t2238 = -0.2e1 * t1921;
t2237 = -0.2e1 * t1922;
t2233 = -g(1) / 0.2e1;
t2232 = g(1) / 0.2e1;
t2231 = -g(2) / 0.2e1;
t2230 = g(2) / 0.2e1;
t2229 = pkin(1) * g(2);
t2225 = t1749 / 0.2e1;
t2224 = t1750 / 0.2e1;
t2223 = t1751 / 0.2e1;
t2222 = t1803 / 0.2e1;
t2221 = t1805 / 0.2e1;
t2220 = t1807 / 0.2e1;
t2219 = t1810 / 0.2e1;
t2218 = t1812 / 0.2e1;
t2217 = t1814 / 0.2e1;
t2216 = t1867 / 0.2e1;
t2215 = t1868 / 0.2e1;
t2214 = t1869 / 0.2e1;
t2213 = t1873 / 0.2e1;
t2212 = t1874 / 0.2e1;
t2211 = t1875 / 0.2e1;
t1904 = t1916 ^ 2;
t2210 = -0.1e1 + 0.2e1 * t1904;
t2209 = 0.4e1 * t1904 - 0.2e1;
t2187 = MDP(4) * t1916;
t2186 = MDP(5) * t1915;
t2185 = qJ(2,1) * t1688;
t2184 = qJ(2,2) * t1687;
t2183 = qJ(2,3) * t1686;
t2179 = t1686 * t1722;
t2178 = t1686 * t1725;
t2177 = t1686 * t1855;
t2176 = t1687 * t1723;
t2175 = t1687 * t1726;
t2174 = t1687 * t1858;
t2173 = t1688 * t1724;
t2172 = t1688 * t1727;
t2171 = t1688 * t1861;
t1728 = t1731 ^ 2;
t2158 = t1728 * t1890;
t1729 = t1732 ^ 2;
t2157 = t1729 * t1893;
t1730 = t1733 ^ 2;
t2156 = t1730 * t1896;
t2152 = t1749 * t1686;
t2151 = t1750 * t1687;
t2150 = t1751 * t1688;
t2136 = t1855 * t1952;
t2133 = t1858 * t1952;
t2130 = t1861 * t1952;
t1911 = t1932 ^ 2;
t2128 = t1911 * t1686;
t1912 = t1934 ^ 2;
t2127 = t1912 * t1687;
t1913 = t1936 ^ 2;
t2126 = t1913 * t1688;
t2123 = t1915 * t1916;
t2122 = t1923 * t1952;
t2121 = t1926 * t1932;
t2120 = t1930 * t1936;
t2119 = t1934 * t1928;
t2117 = g(2) * t2222 + t1809 * t2232;
t2116 = g(2) * t2221 + t1811 * t2232;
t2115 = g(2) * t2220 + t1813 * t2232;
t2114 = g(1) * t2222 + t1809 * t2231;
t2113 = g(1) * t2221 + t1811 * t2231;
t2112 = g(1) * t2220 + t1813 * t2231;
t2108 = t1904 - 0.1e1 / 0.2e1;
t2107 = t1911 - 0.1e1 / 0.2e1;
t2106 = t1912 - 0.1e1 / 0.2e1;
t2105 = t1913 - 0.1e1 / 0.2e1;
t2095 = -0.8e1 * t2123;
t2094 = 0.4e1 * t2123;
t2078 = t1856 * t2158;
t2077 = t1728 * t1857 * t1889;
t2076 = t1859 * t2157;
t2075 = t1729 * t1860 * t1892;
t2074 = t1862 * t2156;
t2073 = t1730 * t1863 * t1895;
t2057 = t1888 * t2118;
t2056 = t1915 * t2124;
t2055 = t1916 * t2124;
t2054 = t1686 * t2121;
t2053 = t1688 * t2120;
t2052 = t1687 * t2119;
t2051 = t1891 * t2118;
t2050 = t1894 * t2118;
t2041 = pkin(2) * t1728 * t2135;
t2040 = pkin(2) * t1729 * t2132;
t2039 = pkin(2) * t1730 * t2129;
t2032 = t2121 * t2123;
t2031 = t2120 * t2123;
t2030 = t2119 * t2123;
t2029 = -t1749 * t2137 / 0.2e1;
t2027 = -t1750 * t2134 / 0.2e1;
t2025 = -t1751 * t2131 / 0.2e1;
t2017 = t2057 * t2155;
t2016 = t2051 * t2154;
t2015 = t2050 * t2153;
t2014 = t2071 * t2228;
t2013 = t2069 * t2227;
t2012 = t2067 * t2226;
t2005 = g(2) * t2219 + t1804 * t2233;
t2004 = g(2) * t2218 + t1806 * t2233;
t2003 = g(2) * t2217 + t1808 * t2233;
t2002 = g(1) * t2219 + t1804 * t2230;
t2001 = g(1) * t2218 + t1806 * t2230;
t2000 = g(1) * t2217 + t1808 * t2230;
t1999 = pkin(2) * t2017;
t1998 = pkin(2) * t2016;
t1997 = pkin(2) * t2015;
t1996 = t1926 * t2017;
t1995 = t1928 * t2016;
t1994 = t1930 * t2015;
t1990 = t1932 * t1996;
t1989 = t1936 * t1994;
t1988 = t1934 * t1995;
t1656 = -t1683 - t2270 * t1888 + (-t2083 + (-t1728 * qJ(2,3) - t2182 / 0.4e1 + t2169 / 0.8e1) * t1889) * t1856 - t1743;
t1977 = t1686 * t1916 + t2017 * t2265;
t1978 = t1686 * t1915 + t2017 * t2264;
t1671 = t1926 * t1978 - t1932 * t1977;
t1674 = t1926 * t1977 + t1932 * t1978;
t1984 = t1671 * MDP(13) + t1674 * MDP(14) + t1656 * MDP(7);
t1657 = -t1684 - t2271 * t1891 + (-t2081 + (-t1729 * qJ(2,2) - t2181 / 0.4e1 + t2168 / 0.8e1) * t1892) * t1859 - t1744;
t1975 = t1687 * t1916 + t2016 * t2265;
t1976 = t1687 * t1915 + t2016 * t2264;
t1672 = t1928 * t1976 - t1934 * t1975;
t1675 = t1928 * t1975 + t1934 * t1976;
t1983 = t1672 * MDP(13) + t1675 * MDP(14) + t1657 * MDP(7);
t1658 = -t1685 - t2272 * t1894 + (-t2079 + (-t1730 * qJ(2,1) - t2180 / 0.4e1 + t2170 / 0.8e1) * t1895) * t1862 - t1745;
t1973 = t1688 * t1916 + t2015 * t2265;
t1974 = t1688 * t1915 + t2015 * t2264;
t1673 = t1930 * t1974 - t1936 * t1973;
t1676 = t1930 * t1973 + t1936 * t1974;
t1982 = t1673 * MDP(13) + t1676 * MDP(14) + t1658 * MDP(7);
t1752 = t1857 * t1867 * t2124 + t1855 * t2122;
t1710 = t1752 * t2239 + t2017 * t2266;
t1641 = t1686 * t2208 + t1864 * t1999 + t2269 * t2213 + t1710 * t2216 - pkin(2) * (-t1686 * t1932 - t1996) - t2005 + t2114;
t1644 = t1870 * t1999 - t2269 * t2216 + t1710 * t2213 - t2002 + t2117 + (-t2111 * t1686 + t1932 * t2017) * pkin(2);
t1941 = pkin(1) * g(1);
t1849 = -g(2) * qJ(2,3) + t1941;
t1850 = g(1) * qJ(2,3) + t2229;
t1647 = (t1849 * t1927 - t1850 * t1933) * t1884 + (t1849 * t1933 + t1850 * t1927) * t1881 + (t2014 + t2183) * qJ(2,3) + (t1683 + t2260) * pkin(1);
t1650 = (-t2128 / 0.2e1 + t1755 / 0.4e1 - t1758 / 0.4e1 - t1767 / 0.4e1 - t1779 / 0.4e1) * t2095 + t2209 * t2054 + ((-t1695 * t1889 / 0.8e1 + t1890 * t2087) * t2095 - ((0.8e1 * t1911 - 0.4e1) * t1904 - 0.8e1 * t2032 - 0.4e1 * t1911 + 0.2e1) * t2057) * t2155;
t1653 = (t1686 - 0.4e1 * t1990 - 0.2e1 * t2128) * t1904 + (t2054 / 0.2e1 - t2107 * t2017) * t2094 + 0.2e1 * t1990 + t2128;
t1680 = t2014 + 0.2e1 * t2183 - t1746;
t1734 = -t1752 * t1916 + t1856 * t2056;
t1737 = t1752 * t1915 + t1856 * t2055;
t1713 = -t1734 * t1926 + t1737 * t1932;
t1716 = -t1734 * t1932 - t1737 * t1926;
t1969 = t1686 * MDP(1) + t1713 * MDP(10) + t1716 * MDP(11) + t1641 * MDP(13) + t1644 * MDP(14) + t1743 * MDP(2) + t1746 * MDP(3) + t1680 * MDP(6) + t1647 * MDP(7) + t1653 * MDP(8) + t1650 * MDP(9);
t1753 = t1860 * t1868 * t2124 + t1858 * t2122;
t1711 = t1753 * t2238 + t2016 * t2266;
t1642 = t1687 * t2207 + t1865 * t1998 + t2268 * t2212 + t1711 * t2215 - pkin(2) * (-t1687 * t1934 - t1995) - t2004 + t2113;
t1645 = t1871 * t1998 - t2268 * t2215 + t1711 * t2212 - t2001 + t2116 + (-t2110 * t1687 + t1934 * t2016) * pkin(2);
t1851 = -g(2) * qJ(2,2) + t1941;
t1852 = g(1) * qJ(2,2) + t2229;
t1648 = (t1851 * t1929 - t1852 * t1935) * t1885 + (t1851 * t1935 + t1852 * t1929) * t1882 + (t2013 + t2184) * qJ(2,2) + (t1684 + t2261) * pkin(1);
t1651 = (-t2127 / 0.2e1 + t1756 / 0.4e1 - t1759 / 0.4e1 - t1768 / 0.4e1 - t1780 / 0.4e1) * t2095 + t2209 * t2052 + ((-t1696 * t1892 / 0.8e1 + t1893 * t2086) * t2095 - ((0.8e1 * t1912 - 0.4e1) * t1904 - 0.8e1 * t2030 - 0.4e1 * t1912 + 0.2e1) * t2051) * t2154;
t1655 = (t1687 - 0.4e1 * t1988 - 0.2e1 * t2127) * t1904 + (t2052 / 0.2e1 - t2106 * t2016) * t2094 + 0.2e1 * t1988 + t2127;
t1681 = t2013 + 0.2e1 * t2184 - t1747;
t1735 = -t1753 * t1916 + t1859 * t2056;
t1738 = t1753 * t1915 + t1859 * t2055;
t1714 = -t1735 * t1928 + t1738 * t1934;
t1717 = -t1735 * t1934 - t1738 * t1928;
t1968 = t1687 * MDP(1) + t1714 * MDP(10) + t1717 * MDP(11) + t1642 * MDP(13) + t1645 * MDP(14) + t1744 * MDP(2) + t1747 * MDP(3) + t1681 * MDP(6) + t1648 * MDP(7) + t1655 * MDP(8) + t1651 * MDP(9);
t1754 = t1863 * t1869 * t2124 + t1861 * t2122;
t1712 = t1754 * t2237 + t2015 * t2266;
t1643 = t1688 * t2206 + t1866 * t1997 + t2267 * t2211 + t1712 * t2214 - pkin(2) * (-t1688 * t1936 - t1994) - t2003 + t2112;
t1646 = t1872 * t1997 - t2267 * t2214 + t1712 * t2211 - t2000 + t2115 + (-t2109 * t1688 + t1936 * t2015) * pkin(2);
t1853 = -g(2) * qJ(2,1) + t1941;
t1854 = g(1) * qJ(2,1) + t2229;
t1649 = (t1853 * t1931 - t1854 * t1937) * t1886 + (t1853 * t1937 + t1854 * t1931) * t1883 + (t2012 + t2185) * qJ(2,1) + (t1685 + t2262) * pkin(1);
t1652 = (-t2126 / 0.2e1 + t1757 / 0.4e1 - t1760 / 0.4e1 - t1769 / 0.4e1 - t1781 / 0.4e1) * t2095 + t2209 * t2053 + ((-t1697 * t1895 / 0.8e1 + t1896 * t2085) * t2095 - ((0.8e1 * t1913 - 0.4e1) * t1904 - 0.8e1 * t2031 - 0.4e1 * t1913 + 0.2e1) * t2050) * t2153;
t1654 = (t1688 - 0.4e1 * t1989 - 0.2e1 * t2126) * t1904 + (t2053 / 0.2e1 - t2105 * t2015) * t2094 + 0.2e1 * t1989 + t2126;
t1682 = t2012 + 0.2e1 * t2185 - t1748;
t1736 = -t1754 * t1916 + t1862 * t2056;
t1739 = t1754 * t1915 + t1862 * t2055;
t1715 = -t1736 * t1930 + t1739 * t1936;
t1718 = -t1736 * t1936 - t1739 * t1930;
t1967 = t1688 * MDP(1) + t1715 * MDP(10) + t1718 * MDP(11) + t1643 * MDP(13) + t1646 * MDP(14) + t1745 * MDP(2) + t1748 * MDP(3) + t1682 * MDP(6) + t1649 * MDP(7) + t1654 * MDP(8) + t1652 * MDP(9);
t1679 = t1688 * t2237 + (-t1709 * t1733 + t1730 * t2250) * t2129;
t1678 = t1687 * t2238 + (-t1708 * t1732 + t1729 * t2250) * t2132;
t1677 = t1686 * t2239 + (-t1707 * t1731 + t1728 * t2250) * t2135;
t1664 = g(1) * t1775 + g(2) * t2006 + 0.2e1 * t1685 + t2262;
t1663 = g(1) * t1774 + g(2) * t2007 + 0.2e1 * t1684 + t2261;
t1662 = g(1) * t1773 + g(2) * t2008 + 0.2e1 * t1683 + t2260;
t1 = [(t1722 * t2078 + t1723 * t2076 + t1724 * t2074) * MDP(6) + (t1925 - g(1)) * MDP(15) - (t1982 * t1724 - t1967 * t2006) * t1894 - (t1983 * t1723 - t1968 * t2007) * t1891 - (t1984 * t1722 - t1969 * t2008) * t1888 + (-(-t2006 * t2247 - t2173) * t1894 - (-t2007 * t2248 - t2176) * t1891 - (-t2008 * t2249 - t2179) * t1888) * t2187 + (-(t1664 * t2006 + t2173) * t1894 - (t1663 * t2007 + t2176) * t1891 - (t1662 * t2008 + t2179) * t1888) * t2186; (t1725 * t2078 + t1726 * t2076 + t1727 * t2074) * MDP(6) + (t1924 - g(2)) * MDP(15) - (t1982 * t1727 + t1967 * t1775) * t1894 - (t1983 * t1726 + t1968 * t1774) * t1891 - (t1984 * t1725 + t1969 * t1773) * t1888 + (-(t1775 * t2247 - t2172) * t1894 - (t1774 * t2248 - t2175) * t1891 - (t1773 * t2249 - t2178) * t1888) * t2187 + (-(-t1664 * t1775 + t2172) * t1894 - (-t1663 * t1774 + t2175) * t1891 - (-t1662 * t1773 + t2178) * t1888) * t2186; (-t1686 * t2063 - t1687 * t2061 - t1688 * t2059) * MDP(1) + (-t1743 * t2063 - t1744 * t2061 - t1745 * t2059) * MDP(2) + (-t1746 * t2063 - t1747 * t2061 - t1748 * t2059) * MDP(3) + (-(t1661 * t2241 - t2150 / 0.2e1) * t2131 - (t1660 * t2242 - t2151 / 0.2e1) * t2134 - (t1659 * t2243 - t2152 / 0.2e1) * t2137) * t2187 + (-(-t1869 * t1664 + t2150 / 0.2e1) * t2131 - (-t1868 * t1663 + t2151 / 0.2e1) * t2134 - (-t1867 * t1662 + t2152 / 0.2e1) * t2137) * t2186 + (t1857 * t2158 * t2225 + t1860 * t2157 * t2224 + t1863 * t2156 * t2223 - t1680 * t2063 - t1681 * t2061 - t1682 * t2059) * MDP(6) + (-(t1869 * t1649 + t1658 * t2223) * t2131 - (t1868 * t1648 + t1657 * t2224) * t2134 - (t1867 * t1647 + t1656 * t2225) * t2137) * MDP(7) + (-t1653 * t2063 - t1654 * t2059 - t1655 * t2061 - 0.2e1 * ((t2105 * t2123 + t2108 * t2120) * t2073 + (t2106 * t2123 + t2108 * t2119) * t2075 + (t2107 * t2123 + t2108 * t2121) * t2077) * t1952) * MDP(8) + (-t1650 * t2063 - t1651 * t2061 - t1652 * t2059 + ((-t1913 * t2209 + 0.4e1 * t2031 + t2210) * t2073 + (-t1912 * t2209 + 0.4e1 * t2030 + t2210) * t2075 + (-t1911 * t2209 + 0.4e1 * t2032 + t2210) * t2077) * t1952) * MDP(9) + (-t1713 * t2063 - t1714 * t2061 - t1715 * t2059 + ((t1915 * t1936 + t1916 * t1930) * t2171 + (t1915 * t1934 + t1916 * t1928) * t2174 + (t1915 * t1932 + t1916 * t1926) * t2177) * t1952) * MDP(10) + (-t1716 * t2063 - t1717 * t2061 - t1718 * t2059 + ((-t1915 * t1930 + t1916 * t1936) * t2171 + (-t1915 * t1928 + t1916 * t1934) * t2174 + (-t1915 * t1926 + t1916 * t1932) * t2177) * t1952) * MDP(11) + (t1752 * t1855 + t1753 * t1858 + t1754 * t1861) * t1952 * MDP(12) + (-t1643 * t2059 + t1673 * t2025 + (t1679 * t2214 - g(3) * t1875 + (t1866 / 0.2e1 + t1930 / 0.2e1) * t2039 + t2003 + t2112) * t2130 - t1642 * t2061 + t1672 * t2027 + (t1678 * t2215 - g(3) * t1874 + (t1865 / 0.2e1 + t1928 / 0.2e1) * t2040 + t2004 + t2113) * t2133 - t1641 * t2063 + t1671 * t2029 + (t1677 * t2216 - g(3) * t1873 + (t1864 / 0.2e1 + t1926 / 0.2e1) * t2041 + t2005 + t2114) * t2136) * MDP(13) + (-t1646 * t2059 + t1676 * t2025 + (t1679 * t2211 + g(3) * t1869 + (t1872 / 0.2e1 + t1936 / 0.2e1) * t2039 + t2000 + t2115) * t2130 - t1645 * t2061 + t1675 * t2027 + (t1678 * t2212 + g(3) * t1868 + (t1871 / 0.2e1 + t1934 / 0.2e1) * t2040 + t2001 + t2116) * t2133 - t1644 * t2063 + t1674 * t2029 + (t1677 * t2213 + g(3) * t1867 + (t1870 / 0.2e1 + t1932 / 0.2e1) * t2041 + t2002 + t2117) * t2136) * MDP(14) + (t1923 - g(3)) * MDP(15);];
tauX  = t1;
