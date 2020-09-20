% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR9V1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR9V1G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:57:30
% EndTime: 2020-08-06 18:57:31
% DurationCPUTime: 1.59s
% Computational Cost: add. (852->175), mult. (1080->310), div. (153->10), fcn. (1101->26), ass. (0->161)
t1923 = pkin(7) + qJ(3,3);
t1897 = sin(t1923);
t1900 = cos(t1923);
t1931 = legFrame(3,2);
t1907 = sin(t1931);
t1910 = cos(t1931);
t1941 = cos(qJ(1,3));
t1997 = t1910 * t1941;
t1852 = t1907 * t1897 + t1900 * t1997;
t1924 = pkin(7) + qJ(3,2);
t1898 = sin(t1924);
t1901 = cos(t1924);
t1932 = legFrame(2,2);
t1908 = sin(t1932);
t1911 = cos(t1932);
t1943 = cos(qJ(1,2));
t1996 = t1911 * t1943;
t1854 = t1908 * t1898 + t1901 * t1996;
t1925 = pkin(7) + qJ(3,1);
t1899 = sin(t1925);
t1902 = cos(t1925);
t1933 = legFrame(1,2);
t1909 = sin(t1933);
t1912 = cos(t1933);
t1945 = cos(qJ(1,1));
t1995 = t1912 * t1945;
t1856 = t1909 * t1899 + t1902 * t1995;
t1896 = 0.1e1 / t1902;
t2034 = pkin(5) + pkin(6);
t1921 = qJ(2,1) + t2034;
t1906 = 0.1e1 / t1921;
t2003 = t1896 * t1906;
t1880 = t1912 * g(1) - t1909 * g(2);
t1939 = sin(qJ(1,1));
t2037 = g(3) * t1945 + t1880 * t1939;
t2038 = t2037 * t2003;
t1879 = t1911 * g(1) - t1908 * g(2);
t1937 = sin(qJ(1,2));
t1868 = g(3) * t1943 + t1879 * t1937;
t1895 = 0.1e1 / t1901;
t1920 = qJ(2,2) + t2034;
t1905 = 0.1e1 / t1920;
t2006 = t1895 * t1905;
t2040 = t1868 * t2006;
t1878 = t1910 * g(1) - t1907 * g(2);
t1935 = sin(qJ(1,3));
t1866 = g(3) * t1941 + t1878 * t1935;
t1894 = 0.1e1 / t1900;
t1919 = qJ(2,3) + t2034;
t1904 = 0.1e1 / t1919;
t2009 = t1894 * t1904;
t2042 = t1866 * t2009;
t2045 = t1852 * t2042 + t1854 * t2040 + t1856 * t2038;
t2000 = t1907 * t1941;
t1851 = t1910 * t1897 - t1900 * t2000;
t1999 = t1908 * t1943;
t1853 = t1911 * t1898 - t1901 * t1999;
t1998 = t1909 * t1945;
t1855 = t1912 * t1899 - t1902 * t1998;
t2044 = t1851 * t2042 + t1853 * t2040 + t1855 * t2038;
t1989 = t1939 * t1906;
t1977 = t2037 * t1989;
t1991 = t1937 * t1905;
t2041 = t1868 * t1991;
t1992 = t1935 * t1904;
t2043 = t1866 * t1992;
t2039 = -t2041 - t2043 - t1977;
t2036 = 2 * pkin(3);
t1930 = cos(pkin(7));
t2035 = 0.2e1 * t1930 ^ 2;
t1940 = cos(qJ(3,3));
t1926 = t1940 ^ 2;
t2030 = t1926 * pkin(3);
t1942 = cos(qJ(3,2));
t1927 = t1942 ^ 2;
t2029 = t1927 * pkin(3);
t1944 = cos(qJ(3,1));
t1928 = t1944 ^ 2;
t2028 = t1928 * pkin(3);
t2027 = t1940 * pkin(2);
t2026 = t1942 * pkin(2);
t2025 = t1944 * pkin(2);
t1917 = pkin(1) * t1941;
t1842 = t1878 * (t1935 * pkin(1) - t1941 * qJ(2,3)) + g(3) * (t1935 * qJ(2,3) + t1917);
t2024 = t1842 * t1894;
t1918 = pkin(1) * t1943;
t1843 = t1879 * (t1937 * pkin(1) - t1943 * qJ(2,2)) + g(3) * (t1937 * qJ(2,2) + t1918);
t2023 = t1843 * t1895;
t1913 = t1945 * pkin(1);
t1844 = t1880 * (t1939 * pkin(1) - t1945 * qJ(2,1)) + g(3) * (t1939 * qJ(2,1) + t1913);
t2022 = t1844 * t1896;
t2021 = t2037 * t1906;
t1929 = sin(pkin(7));
t1934 = sin(qJ(3,3));
t1994 = t1929 * t1934;
t2020 = t1866 / (t1930 * t1940 - t1994);
t2019 = t1866 * t1904;
t1936 = sin(qJ(3,2));
t1993 = t1929 * t1936;
t2018 = t1868 / (t1930 * t1942 - t1993);
t2017 = t1868 * t1905;
t1938 = sin(qJ(3,1));
t1990 = t1938 * t1929;
t2016 = t2037 / (t1944 * t1930 - t1990);
t1947 = pkin(2) / 0.2e1;
t2012 = (t1940 * pkin(3) + t1947) * t1934;
t2011 = (t1942 * pkin(3) + t1947) * t1936;
t2010 = (t1944 * pkin(3) + t1947) * t1938;
t2008 = t1894 * t1907;
t2007 = t1894 * t1910;
t2005 = t1895 * t1908;
t2004 = t1895 * t1911;
t2002 = t1896 * t1909;
t2001 = t1896 * t1912;
t1903 = pkin(1) * t1929;
t1988 = t1940 * (-t1934 * pkin(3) + t1903);
t1987 = t1942 * (-t1936 * pkin(3) + t1903);
t1986 = t1944 * (-t1938 * pkin(3) + t1903);
t1985 = t1939 * t1921 + t1913;
t1984 = t1935 * t1919 + t1917;
t1983 = t1937 * t1920 + t1918;
t1965 = -g(3) * t1935 + t1878 * t1941;
t1975 = t1965 * t2009;
t1964 = -g(3) * t1937 + t1879 * t1943;
t1973 = t1964 * t2006;
t1963 = -g(3) * t1939 + t1880 * t1945;
t1971 = t1963 * t2003;
t1970 = t1941 * t1994;
t1969 = t1943 * t1993;
t1968 = t1945 * t1990;
t1960 = t1897 * t2042;
t1959 = t1898 * t2040;
t1958 = t1899 * t2038;
t1957 = pkin(2) * t1970 + (t1970 * t2036 - t1984) * t1940;
t1956 = pkin(2) * t1969 + (t1969 * t2036 - t1983) * t1942;
t1955 = pkin(2) * t1968 + (t1968 * t2036 - t1985) * t1944;
t1951 = t1963 * t1989 + t1964 * t1991 + t1965 * t1992;
t1950 = t1851 * t1975 + t1853 * t1973 + t1855 * t1971;
t1949 = t1852 * t1975 + t1854 * t1973 + t1856 * t1971;
t1948 = 1 / pkin(3);
t1946 = -pkin(3) / 0.2e1;
t1893 = t1930 * pkin(2) + pkin(1);
t1883 = t2028 + t2025 / 0.2e1 + t1946;
t1882 = t2029 + t2026 / 0.2e1 + t1946;
t1881 = t2030 + t2027 / 0.2e1 + t1946;
t1877 = t1909 * g(1) + t1912 * g(2);
t1876 = t1908 * g(1) + t1911 * g(2);
t1875 = t1907 * g(1) + t1910 * g(2);
t1850 = pkin(1) * t1938 + (-pkin(3) + t2025 + 0.2e1 * t2028) * t1929;
t1849 = pkin(1) * t1936 + (-pkin(3) + t2026 + 0.2e1 * t2029) * t1929;
t1848 = pkin(1) * t1934 + (-pkin(3) + t2027 + 0.2e1 * t2030) * t1929;
t1847 = t1985 * t1990 + (t1928 - 0.1e1) * t1945 * pkin(3);
t1846 = t1983 * t1993 + (t1927 - 0.1e1) * t1943 * pkin(3);
t1845 = t1984 * t1994 + (t1926 - 0.1e1) * t1941 * pkin(3);
t1841 = t1877 * t1899 + t1902 * t1963;
t1840 = -t1877 * t1902 + t1899 * t1963;
t1839 = t1876 * t1898 + t1901 * t1964;
t1838 = -t1876 * t1901 + t1898 * t1964;
t1837 = t1875 * t1897 + t1900 * t1965;
t1836 = -t1875 * t1900 + t1897 * t1965;
t1 = [0, t2045, t1949, t2045 * t1930, -t2045 * t1929, -t1949, (t1856 * t2022 - ((t1883 * t1995 + t1909 * t2010) * t2035 + (t1909 * t1850 - t1955 * t1912) * t1930 - t1847 * t1912 + t1909 * t1986) * t2016) * t1906 + (t1854 * t2023 - ((t1882 * t1996 + t1908 * t2011) * t2035 + (t1908 * t1849 - t1956 * t1911) * t1930 - t1846 * t1911 + t1908 * t1987) * t2018) * t1905 + (t1852 * t2024 - ((t1881 * t1997 + t1907 * t2012) * t2035 + (t1907 * t1848 - t1957 * t1910) * t1930 - t1845 * t1910 + t1907 * t1988) * t2020) * t1904, 0, 0, 0, 0, 0, t1852 * t2019 + t1854 * t2017 + t1856 * t2021 + (t1836 * t2008 + t1838 * t2005 + t1840 * t2002) * t1948, -t1852 * t1960 - t1854 * t1959 - t1856 * t1958 + (t1837 * t2008 + t1839 * t2005 + t1841 * t2002) * t1948, -g(1); 0, t2044, t1950, t2044 * t1930, -t2044 * t1929, -t1950, (t1855 * t2022 - ((-t1883 * t1998 + t1912 * t2010) * t2035 + (t1912 * t1850 + t1955 * t1909) * t1930 + t1847 * t1909 + t1912 * t1986) * t2016) * t1906 + (t1853 * t2023 - ((-t1882 * t1999 + t1911 * t2011) * t2035 + (t1911 * t1849 + t1956 * t1908) * t1930 + t1846 * t1908 + t1911 * t1987) * t2018) * t1905 + (t1851 * t2024 - ((-t1881 * t2000 + t1910 * t2012) * t2035 + (t1910 * t1848 + t1957 * t1907) * t1930 + t1845 * t1907 + t1910 * t1988) * t2020) * t1904, 0, 0, 0, 0, 0, t1851 * t2019 + t1853 * t2017 + t1855 * t2021 + (t1836 * t2007 + t1838 * t2004 + t1840 * t2001) * t1948, -t1851 * t1960 - t1853 * t1959 - t1855 * t1958 + (t1837 * t2007 + t1839 * t2004 + t1841 * t2001) * t1948, -g(2); 0, t2039, -t1951, t2039 * t1930, -t2039 * t1929, t1951, (-t1921 * t1945 * t2037 + (-t1844 - (-pkin(3) * t1902 - t1893) * t2037) * t1939) * t1906 + (-t1920 * t1943 * t1868 + (-t1843 - (-pkin(3) * t1901 - t1893) * t1868) * t1937) * t1905 + (-t1919 * t1941 * t1866 + (-t1842 - (-pkin(3) * t1900 - t1893) * t1866) * t1935) * t1904, 0, 0, 0, 0, 0, -t1900 * t2043 - t1901 * t2041 - t1902 * t1977, t1897 * t2043 + t1898 * t2041 + t1899 * t1977, -g(3);];
tau_reg  = t1;
