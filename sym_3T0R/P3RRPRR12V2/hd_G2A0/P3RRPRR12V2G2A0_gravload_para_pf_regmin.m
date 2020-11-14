% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V2G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:22:32
% EndTime: 2020-08-06 19:22:34
% DurationCPUTime: 2.43s
% Computational Cost: add. (1482->221), mult. (2157->409), div. (180->6), fcn. (1953->18), ass. (0->178)
t2150 = legFrame(3,2);
t2126 = sin(t2150);
t2129 = cos(t2150);
t2108 = t2129 * g(1) - t2126 * g(2);
t2154 = sin(qJ(1,3));
t2160 = cos(qJ(1,3));
t2084 = -g(3) * t2154 + t2108 * t2160;
t2305 = t2084 * t2160;
t2151 = legFrame(2,2);
t2127 = sin(t2151);
t2130 = cos(t2151);
t2109 = t2130 * g(1) - t2127 * g(2);
t2156 = sin(qJ(1,2));
t2162 = cos(qJ(1,2));
t2085 = -g(3) * t2156 + t2109 * t2162;
t2304 = t2085 * t2162;
t2152 = legFrame(1,2);
t2128 = sin(t2152);
t2131 = cos(t2152);
t2110 = t2131 * g(1) - t2128 * g(2);
t2158 = sin(qJ(1,1));
t2164 = cos(qJ(1,1));
t2086 = -g(3) * t2158 + t2110 * t2164;
t2303 = t2086 * t2164;
t2165 = pkin(5) - pkin(6);
t2178 = pkin(1) * t2154 - t2165 * t2160;
t2153 = sin(qJ(2,3));
t2229 = t2153 * t2154;
t2202 = qJ(3,3) * t2229;
t2094 = t2178 + t2202;
t2117 = pkin(1) * t2153 + qJ(3,3);
t2159 = cos(qJ(2,3));
t2147 = t2159 ^ 2;
t2166 = pkin(2) + pkin(3);
t2226 = t2154 * t2166;
t2227 = t2153 * t2166;
t2259 = t2126 * qJ(3,3);
t2052 = (t2129 * t2226 - t2259) * t2147 + (t2094 * t2129 + t2126 * t2227) * t2159 + t2126 * t2117;
t2177 = pkin(1) * t2156 - t2165 * t2162;
t2155 = sin(qJ(2,2));
t2225 = t2155 * t2156;
t2201 = qJ(3,2) * t2225;
t2096 = t2177 + t2201;
t2118 = pkin(1) * t2155 + qJ(3,2);
t2161 = cos(qJ(2,2));
t2148 = t2161 ^ 2;
t2222 = t2156 * t2166;
t2223 = t2155 * t2166;
t2258 = t2127 * qJ(3,2);
t2054 = (t2130 * t2222 - t2258) * t2148 + (t2096 * t2130 + t2127 * t2223) * t2161 + t2127 * t2118;
t2176 = pkin(1) * t2158 - t2165 * t2164;
t2157 = sin(qJ(2,1));
t2221 = t2157 * t2158;
t2200 = qJ(3,1) * t2221;
t2098 = t2176 + t2200;
t2119 = pkin(1) * t2157 + qJ(3,1);
t2163 = cos(qJ(2,1));
t2149 = t2163 ^ 2;
t2218 = t2158 * t2166;
t2219 = t2157 * t2166;
t2257 = t2128 * qJ(3,1);
t2056 = (t2131 * t2218 - t2257) * t2149 + (t2098 * t2131 + t2128 * t2219) * t2163 + t2128 * t2119;
t2211 = t2166 * t2159;
t2123 = t2153 * qJ(3,3);
t2262 = t2123 + pkin(1);
t2175 = t2211 + t2262;
t2102 = 0.1e1 / t2175;
t2210 = t2166 * t2161;
t2124 = t2155 * qJ(3,2);
t2261 = t2124 + pkin(1);
t2174 = t2210 + t2261;
t2103 = 0.1e1 / t2174;
t2209 = t2166 * t2163;
t2125 = t2157 * qJ(3,1);
t2260 = t2125 + pkin(1);
t2173 = t2209 + t2260;
t2104 = 0.1e1 / t2173;
t2107 = t2128 * g(1) + t2131 * g(2);
t2281 = g(3) * t2164 + t2110 * t2158;
t2072 = t2107 * t2157 + t2163 * t2281;
t2169 = 0.1e1 / qJ(3,1);
t2288 = t2072 * t2169;
t2106 = t2127 * g(1) + t2130 * g(2);
t2280 = g(3) * t2162 + t2109 * t2156;
t2068 = t2106 * t2155 + t2161 * t2280;
t2168 = 0.1e1 / qJ(3,2);
t2289 = t2068 * t2168;
t2105 = t2126 * g(1) + t2129 * g(2);
t2279 = g(3) * t2160 + t2108 * t2154;
t2064 = t2105 * t2153 + t2159 * t2279;
t2167 = 0.1e1 / qJ(3,3);
t2290 = t2064 * t2167;
t2298 = t2157 * t2303;
t2299 = t2155 * t2304;
t2300 = t2153 * t2305;
t2302 = (t2052 * t2290 + t2129 * t2300) * t2102 + (t2054 * t2289 + t2130 * t2299) * t2103 + (t2056 * t2288 + t2131 * t2298) * t2104;
t2256 = t2129 * qJ(3,3);
t2051 = (-t2126 * t2226 - t2256) * t2147 + (-t2094 * t2126 + t2129 * t2227) * t2159 + t2129 * t2117;
t2255 = t2130 * qJ(3,2);
t2053 = (-t2127 * t2222 - t2255) * t2148 + (-t2096 * t2127 + t2130 * t2223) * t2161 + t2130 * t2118;
t2254 = t2131 * qJ(3,1);
t2055 = (-t2128 * t2218 - t2254) * t2149 + (-t2098 * t2128 + t2131 * t2219) * t2163 + t2131 * t2119;
t2301 = (t2051 * t2290 - t2126 * t2300) * t2102 + (t2053 * t2289 - t2127 * t2299) * t2103 + (t2055 * t2288 - t2128 * t2298) * t2104;
t2122 = t2158 * t2165;
t2077 = t2173 * t2164 + t2122;
t2242 = t2077 * t2169;
t2121 = t2156 * t2165;
t2076 = t2174 * t2162 + t2121;
t2243 = t2076 * t2168;
t2120 = t2154 * t2165;
t2075 = t2175 * t2160 + t2120;
t2244 = t2075 * t2167;
t2291 = (t2064 * t2159 * t2244 - t2084 * t2229) * t2102 + (t2068 * t2161 * t2243 - t2085 * t2225) * t2103 + (t2072 * t2163 * t2242 - t2086 * t2221) * t2104;
t2272 = -0.2e1 * t2166;
t2063 = -t2105 * t2159 + t2153 * t2279;
t2253 = t2063 * t2167;
t2067 = -t2106 * t2161 + t2155 * t2280;
t2250 = t2067 * t2168;
t2071 = -t2107 * t2163 + t2157 * t2281;
t2247 = t2071 * t2169;
t2241 = (t2178 + 0.2e1 * t2202) * t2166;
t2240 = (t2177 + 0.2e1 * t2201) * t2166;
t2239 = (t2176 + 0.2e1 * t2200) * t2166;
t2238 = t2102 * t2154;
t2237 = t2102 * t2160;
t2236 = t2103 * t2156;
t2235 = t2103 * t2162;
t2234 = t2104 * t2158;
t2233 = t2104 * t2164;
t2232 = (qJ(3,3) + t2166) * (-qJ(3,3) + t2166);
t2231 = (qJ(3,2) + t2166) * (-qJ(3,2) + t2166);
t2230 = (qJ(3,1) + t2166) * (-qJ(3,1) + t2166);
t2132 = t2159 * pkin(2);
t2208 = t2132 + t2123;
t2062 = g(3) * (-t2160 * pkin(5) + (pkin(1) + t2208) * t2154) - t2108 * (t2154 * pkin(5) + (t2132 + t2262) * t2160);
t2217 = t2160 * t2062;
t2133 = t2161 * pkin(2);
t2207 = t2133 + t2124;
t2060 = (-t2162 * pkin(5) + (pkin(1) + t2207) * t2156) * g(3) - t2109 * (t2156 * pkin(5) + (t2133 + t2261) * t2162);
t2216 = t2162 * t2060;
t2134 = t2163 * pkin(2);
t2206 = t2134 + t2125;
t2061 = (-t2164 * pkin(5) + (pkin(1) + t2206) * t2158) * g(3) - (t2158 * pkin(5) + (t2134 + t2260) * t2164) * t2110;
t2215 = t2164 * t2061;
t2214 = t2166 * t2117;
t2213 = t2166 * t2118;
t2212 = t2166 * t2119;
t2205 = qJ(3,1) * t2272;
t2204 = qJ(3,2) * t2272;
t2203 = qJ(3,3) * t2272;
t2196 = t2159 * t2305;
t2195 = t2161 * t2304;
t2194 = t2163 * t2303;
t2193 = t2126 * t2237;
t2192 = t2129 * t2237;
t2191 = t2127 * t2235;
t2190 = t2130 * t2235;
t2189 = t2128 * t2233;
t2188 = t2131 * t2233;
t2181 = t2154 * t2232;
t2180 = t2156 * t2231;
t2179 = t2158 * t2230;
t2172 = t2234 * t2281 + t2236 * t2280 + t2238 * t2279;
t2171 = t2189 * t2281 + t2191 * t2280 + t2193 * t2279;
t2170 = t2188 * t2281 + t2190 * t2280 + t2192 * t2279;
t2101 = pkin(1) * qJ(3,1) - t2157 * t2230;
t2100 = pkin(1) * qJ(3,2) - t2155 * t2231;
t2099 = pkin(1) * qJ(3,3) - t2153 * t2232;
t2092 = t2158 * qJ(3,1) + t2176 * t2157;
t2091 = t2156 * qJ(3,2) + t2177 * t2155;
t2090 = t2154 * qJ(3,3) + t2178 * t2153;
t2059 = -t2107 * t2206 + t2281 * (t2157 * pkin(2) - t2163 * qJ(3,1));
t2058 = -t2106 * t2207 + t2280 * (t2155 * pkin(2) - t2161 * qJ(3,2));
t2057 = -t2105 * t2208 + t2279 * (t2153 * pkin(2) - t2159 * qJ(3,3));
t2050 = (t2071 * t2242 + t2086 * t2158) * t2163 * t2104 + (t2067 * t2243 + t2085 * t2156) * t2161 * t2103 + (t2063 * t2244 + t2084 * t2154) * t2159 * t2102;
t2049 = (t2056 * t2247 - t2131 * t2194) * t2104 + (t2054 * t2250 - t2130 * t2195) * t2103 + (t2052 * t2253 - t2129 * t2196) * t2102;
t2048 = (t2055 * t2247 + t2128 * t2194) * t2104 + (t2053 * t2250 + t2127 * t2195) * t2103 + (t2051 * t2253 + t2126 * t2196) * t2102;
t1 = [0, -t2084 * t2192 - t2085 * t2190 - t2086 * t2188, t2170, 0, 0, 0, 0, 0, t2049, t2302, t2049, -t2170, -t2302, (t2131 * t2215 + (t2056 * t2059 - ((t2128 * t2205 + t2131 * t2179) * t2149 + (-t2128 * t2101 + t2131 * t2239) * t2163 + t2092 * t2254 + t2128 * t2212) * t2071) * t2169) * t2104 + (t2130 * t2216 + (t2054 * t2058 - ((t2127 * t2204 + t2130 * t2180) * t2148 + (-t2127 * t2100 + t2130 * t2240) * t2161 + t2091 * t2255 + t2127 * t2213) * t2067) * t2168) * t2103 + (t2129 * t2217 + (t2052 * t2057 - ((t2126 * t2203 + t2129 * t2181) * t2147 + (-t2126 * t2099 + t2129 * t2241) * t2159 + t2090 * t2256 + t2126 * t2214) * t2063) * t2167) * t2102, -g(1); 0, t2084 * t2193 + t2085 * t2191 + t2086 * t2189, -t2171, 0, 0, 0, 0, 0, t2048, t2301, t2048, t2171, -t2301, (-t2128 * t2215 + (t2055 * t2059 - ((-t2128 * t2179 + t2131 * t2205) * t2149 + (-t2131 * t2101 - t2128 * t2239) * t2163 - t2092 * t2257 + t2131 * t2212) * t2071) * t2169) * t2104 + (-t2127 * t2216 + (t2053 * t2058 - ((-t2127 * t2180 + t2130 * t2204) * t2148 + (-t2130 * t2100 - t2127 * t2240) * t2161 - t2091 * t2258 + t2130 * t2213) * t2067) * t2168) * t2103 + (-t2126 * t2217 + (t2051 * t2057 - ((-t2126 * t2181 + t2129 * t2203) * t2147 + (-t2129 * t2099 - t2126 * t2241) * t2159 - t2090 * t2259 + t2129 * t2214) * t2063) * t2167) * t2102, -g(2); 0, t2084 * t2238 + t2085 * t2236 + t2086 * t2234, -t2172, 0, 0, 0, 0, 0, t2050, t2291, t2050, t2172, -t2291, (-t2158 * t2061 + (t2077 * t2163 * t2059 - (t2164 * t2149 * t2230 + ((0.2e1 * t2125 + pkin(1)) * t2164 + t2122) * t2209 + qJ(3,1) * (t2119 * t2164 + t2122 * t2157)) * t2071) * t2169) * t2104 + (-t2156 * t2060 + (t2076 * t2161 * t2058 - (t2162 * t2148 * t2231 + ((0.2e1 * t2124 + pkin(1)) * t2162 + t2121) * t2210 + qJ(3,2) * (t2118 * t2162 + t2121 * t2155)) * t2067) * t2168) * t2103 + (-t2154 * t2062 + (t2075 * t2159 * t2057 - (t2160 * t2147 * t2232 + ((0.2e1 * t2123 + pkin(1)) * t2160 + t2120) * t2211 + qJ(3,3) * (t2117 * t2160 + t2120 * t2153)) * t2063) * t2167) * t2102, -g(3);];
tau_reg  = t1;
