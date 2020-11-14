% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V2G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:16:19
% EndTime: 2020-08-06 19:16:21
% DurationCPUTime: 2.31s
% Computational Cost: add. (1245->192), mult. (1791->330), div. (138->6), fcn. (1761->18), ass. (0->152)
t2009 = legFrame(3,3);
t1991 = sin(t2009);
t1994 = cos(t2009);
t2013 = sin(qJ(1,3));
t2019 = cos(qJ(1,3));
t1935 = t1991 * t2019 + t1994 * t2013;
t2012 = sin(qJ(2,3));
t1988 = t2012 * qJ(3,3);
t1971 = t1988 + pkin(1);
t2024 = pkin(5) - pkin(6);
t1985 = t2024 * t2019;
t1941 = t1971 * t2013 - t1985;
t1982 = t2013 * t2024;
t2036 = t1971 * t2019 + t1982;
t2018 = cos(qJ(2,3));
t2025 = pkin(2) + pkin(3);
t2061 = t2018 * t2025;
t1893 = t1935 * t2061 + t1941 * t1994 + t1991 * t2036;
t2010 = legFrame(2,3);
t1992 = sin(t2010);
t1995 = cos(t2010);
t2015 = sin(qJ(1,2));
t2021 = cos(qJ(1,2));
t1937 = t1992 * t2021 + t1995 * t2015;
t2014 = sin(qJ(2,2));
t1989 = t2014 * qJ(3,2);
t1973 = t1989 + pkin(1);
t1986 = t2024 * t2021;
t1943 = t1973 * t2015 - t1986;
t1983 = t2015 * t2024;
t2034 = t1973 * t2021 + t1983;
t2020 = cos(qJ(2,2));
t2059 = t2020 * t2025;
t1895 = t1937 * t2059 + t1943 * t1995 + t1992 * t2034;
t2011 = legFrame(1,3);
t1993 = sin(t2011);
t1996 = cos(t2011);
t2017 = sin(qJ(1,1));
t2023 = cos(qJ(1,1));
t1939 = t1993 * t2023 + t1996 * t2017;
t2016 = sin(qJ(2,1));
t1990 = t2016 * qJ(3,1);
t1975 = t1990 + pkin(1);
t1987 = t2024 * t2023;
t1945 = t1975 * t2017 - t1987;
t1984 = t2017 * t2024;
t2032 = t1975 * t2023 + t1984;
t2022 = cos(qJ(2,1));
t2057 = t2022 * t2025;
t1897 = t1939 * t2057 + t1945 * t1996 + t1993 * t2032;
t1934 = t1991 * t2013 - t1994 * t2019;
t1936 = t1992 * t2015 - t1995 * t2021;
t1938 = t1993 * t2017 - t1996 * t2023;
t1952 = 0.1e1 / (t2061 + t1971);
t1953 = 0.1e1 / (t2059 + t1973);
t1954 = 0.1e1 / (t2057 + t1975);
t1948 = g(1) * t1993 - t1996 * g(2);
t1951 = g(1) * t1996 + g(2) * t1993;
t1921 = t1948 * t2023 + t1951 * t2017;
t2074 = t1921 * t2016;
t1947 = g(1) * t1992 - t1995 * g(2);
t1950 = g(1) * t1995 + g(2) * t1992;
t1920 = t1947 * t2021 + t1950 * t2015;
t2075 = t1920 * t2014;
t1946 = g(1) * t1991 - t1994 * g(2);
t1949 = g(1) * t1994 + g(2) * t1991;
t1919 = t1946 * t2019 + t1949 * t2013;
t2076 = t1919 * t2012;
t1930 = (g(1) * t2017 - g(2) * t2023) * t1993 - (g(1) * t2023 + g(2) * t2017) * t1996;
t1915 = -g(3) * t2016 + t1930 * t2022;
t2028 = 0.1e1 / qJ(3,1);
t2095 = t1915 * t2022 * t2028;
t1928 = (g(1) * t2015 - g(2) * t2021) * t1992 - (g(1) * t2021 + g(2) * t2015) * t1995;
t1911 = -g(3) * t2014 + t1928 * t2020;
t2027 = 0.1e1 / qJ(3,2);
t2096 = t1911 * t2020 * t2027;
t1926 = (g(1) * t2013 - g(2) * t2019) * t1991 - (g(1) * t2019 + g(2) * t2013) * t1994;
t1907 = -g(3) * t2012 + t1926 * t2018;
t2026 = 0.1e1 / qJ(3,3);
t2097 = t1907 * t2018 * t2026;
t2105 = (t1893 * t2097 - t1934 * t2076) * t1952 + (t1895 * t2096 - t1936 * t2075) * t1953 + (t1897 * t2095 - t1938 * t2074) * t1954;
t1892 = t1934 * t2061 + t1941 * t1991 - t2036 * t1994;
t1894 = t1936 * t2059 + t1943 * t1992 - t2034 * t1995;
t1896 = t1938 * t2057 + t1945 * t1993 - t2032 * t1996;
t2104 = (t1892 * t2097 + t1935 * t2076) * t1952 + (t1894 * t2096 + t1937 * t2075) * t1953 + (t1896 * t2095 + t1939 * t2074) * t1954;
t2062 = t2016 * t2028;
t2063 = t2014 * t2027;
t2064 = t2012 * t2026;
t2094 = t1907 * t2064 + t1911 * t2063 + t1915 * t2062;
t1906 = g(3) * t2018 + t1926 * t2012;
t1910 = g(3) * t2020 + t1928 * t2014;
t1914 = g(3) * t2022 + t1930 * t2016;
t1916 = t1946 * t2013 - t1949 * t2019;
t1917 = t1947 * t2015 - t1950 * t2021;
t1918 = t1948 * t2017 - t1951 * t2023;
t2068 = t1939 * t1954;
t2070 = t1937 * t1953;
t2072 = t1935 * t1952;
t2093 = t1916 * t2072 + t1917 * t2070 + t1918 * t2068;
t2069 = t1938 * t1954;
t2071 = t1936 * t1953;
t2073 = t1934 * t1952;
t2092 = t1916 * t2073 + t1917 * t2071 + t1918 * t2069;
t2085 = qJ(3,1) * t2022;
t2084 = qJ(3,2) * t2020;
t2083 = qJ(3,3) * t2018;
t2055 = t2018 * pkin(2) + t1988;
t1901 = -g(3) * t2055 - t1916 * (pkin(2) * t2012 - t2083);
t2082 = t1901 * t2018;
t2054 = t2020 * pkin(2) + t1989;
t1902 = -g(3) * t2054 - t1917 * (pkin(2) * t2014 - t2084);
t2081 = t1902 * t2020;
t2053 = t2022 * pkin(2) + t1990;
t1903 = -g(3) * t2053 - t1918 * (pkin(2) * t2016 - t2085);
t2080 = t1903 * t2022;
t2079 = t1906 * t2026;
t2078 = t1910 * t2027;
t2077 = t1914 * t2028;
t2067 = t1952 * t2018;
t2066 = t1953 * t2020;
t2065 = t1954 * t2022;
t2040 = (qJ(3,3) + t2025) * (-qJ(3,3) + t2025) * t2018 ^ 2;
t2039 = (qJ(3,2) + t2025) * (-qJ(3,2) + t2025) * t2020 ^ 2;
t2038 = (qJ(3,1) + t2025) * (-qJ(3,1) + t2025) * t2022 ^ 2;
t1970 = 0.2e1 * t1988 + pkin(1);
t2037 = t1970 * t2019 + t1982;
t1972 = 0.2e1 * t1989 + pkin(1);
t2035 = t1972 * t2021 + t1983;
t1974 = 0.2e1 * t1990 + pkin(1);
t2033 = t1974 * t2023 + t1984;
t1976 = pkin(1) * t2012 + qJ(3,3);
t2031 = t1976 * t2019 + t2012 * t1982;
t1977 = pkin(1) * t2014 + qJ(3,2);
t2030 = t1977 * t2021 + t2014 * t1983;
t1978 = pkin(1) * t2016 + qJ(3,1);
t2029 = t1978 * t2023 + t2016 * t1984;
t1957 = pkin(1) + t2053;
t1956 = pkin(1) + t2054;
t1955 = pkin(1) + t2055;
t1944 = t1974 * t2017 - t1987;
t1942 = t1972 * t2015 - t1986;
t1940 = t1970 * t2013 - t1985;
t1933 = t1978 * t2017 - t2016 * t1987;
t1932 = t1977 * t2015 - t2014 * t1986;
t1931 = t1976 * t2013 - t2012 * t1985;
t1900 = t1951 * (-pkin(5) * t2023 + t1957 * t2017) + t1948 * (pkin(5) * t2017 + t1957 * t2023);
t1899 = t1950 * (-pkin(5) * t2021 + t1956 * t2015) + t1947 * (pkin(5) * t2015 + t1956 * t2021);
t1898 = t1949 * (-pkin(5) * t2019 + t1955 * t2013) + t1946 * (pkin(5) * t2013 + t1955 * t2019);
t1891 = -t1906 * t2064 - t1910 * t2063 - t1914 * t2062;
t1890 = (t1896 * t2077 - t1921 * t1939) * t2065 + (t1894 * t2078 - t1920 * t1937) * t2066 + (t1892 * t2079 - t1919 * t1935) * t2067;
t1889 = (-t1897 * t2077 - t1921 * t1938) * t2065 + (-t1895 * t2078 - t1920 * t1936) * t2066 + (-t1893 * t2079 - t1919 * t1934) * t2067;
t1 = [0, -t1919 * t2072 - t1920 * t2070 - t1921 * t2068, t2093, 0, 0, 0, 0, 0, t1890, t2104, t1890, -t2093, -t2104, (-t1939 * t1900 + (-t1896 * t2080 + (-t1938 * t2038 - (t1993 * t1944 - t2033 * t1996) * t2057 - (t1993 * t1933 - t2029 * t1996) * qJ(3,1)) * t1914) * t2028) * t1954 + (-t1937 * t1899 + (-t1894 * t2081 + (-t1936 * t2039 - (t1992 * t1942 - t2035 * t1995) * t2059 - (t1992 * t1932 - t2030 * t1995) * qJ(3,2)) * t1910) * t2027) * t1953 + (-t1935 * t1898 + (-t1892 * t2082 + (-t1934 * t2040 - (t1991 * t1940 - t2037 * t1994) * t2061 - (t1991 * t1931 - t2031 * t1994) * qJ(3,3)) * t1906) * t2026) * t1952, -g(1); 0, -t1919 * t2073 - t1920 * t2071 - t1921 * t2069, t2092, 0, 0, 0, 0, 0, t1889, -t2105, t1889, -t2092, t2105, (-t1938 * t1900 + (t1897 * t2080 + (t1939 * t2038 + (t1944 * t1996 + t2033 * t1993) * t2057 + (t1933 * t1996 + t2029 * t1993) * qJ(3,1)) * t1914) * t2028) * t1954 + (-t1936 * t1899 + (t1895 * t2081 + (t1937 * t2039 + (t1942 * t1995 + t2035 * t1992) * t2059 + (t1932 * t1995 + t2030 * t1992) * qJ(3,2)) * t1910) * t2027) * t1953 + (-t1934 * t1898 + (t1893 * t2082 + (t1935 * t2040 + (t1940 * t1994 + t2037 * t1991) * t2061 + (t1931 * t1994 + t2031 * t1991) * qJ(3,3)) * t1906) * t2026) * t1952, -g(2); 0, 0, 0, 0, 0, 0, 0, 0, t1891, -t2094, t1891, 0, t2094, (t2016 * t1903 + (t2016 * t2025 - t2085) * t1914) * t2028 + (t2014 * t1902 + (t2014 * t2025 - t2084) * t1910) * t2027 + (t2012 * t1901 + (t2012 * t2025 - t2083) * t1906) * t2026, -g(3);];
tau_reg  = t1;
