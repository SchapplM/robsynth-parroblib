% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V2G3A0
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
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V2G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:29:06
% EndTime: 2020-08-06 19:29:09
% DurationCPUTime: 2.53s
% Computational Cost: add. (1482->215), mult. (2157->409), div. (180->6), fcn. (1953->18), ass. (0->175)
t2171 = sin(qJ(2,3));
t2178 = cos(qJ(1,3));
t2256 = t2171 * t2178;
t2224 = qJ(3,3) * t2256;
t2172 = sin(qJ(1,3));
t2183 = pkin(5) - pkin(6);
t2232 = pkin(1) * t2178 + t2172 * t2183;
t2106 = t2224 + t2232;
t2132 = pkin(1) * t2171 + qJ(3,3);
t2168 = legFrame(3,2);
t2141 = sin(t2168);
t2144 = cos(t2168);
t2177 = cos(qJ(2,3));
t2165 = t2177 ^ 2;
t2184 = pkin(2) + pkin(3);
t2247 = t2178 * t2184;
t2255 = t2171 * t2184;
t2292 = t2141 * qJ(3,3);
t2064 = (t2144 * t2247 - t2292) * t2165 + (t2106 * t2144 + t2141 * t2255) * t2177 + t2141 * t2132;
t2173 = sin(qJ(2,2));
t2180 = cos(qJ(1,2));
t2253 = t2173 * t2180;
t2225 = qJ(3,2) * t2253;
t2174 = sin(qJ(1,2));
t2231 = pkin(1) * t2180 + t2174 * t2183;
t2108 = t2225 + t2231;
t2133 = pkin(1) * t2173 + qJ(3,2);
t2169 = legFrame(2,2);
t2142 = sin(t2169);
t2145 = cos(t2169);
t2179 = cos(qJ(2,2));
t2166 = t2179 ^ 2;
t2246 = t2180 * t2184;
t2252 = t2173 * t2184;
t2291 = t2142 * qJ(3,2);
t2066 = (t2145 * t2246 - t2291) * t2166 + (t2108 * t2145 + t2142 * t2252) * t2179 + t2142 * t2133;
t2175 = sin(qJ(2,1));
t2182 = cos(qJ(1,1));
t2250 = t2175 * t2182;
t2226 = qJ(3,1) * t2250;
t2176 = sin(qJ(1,1));
t2230 = pkin(1) * t2182 + t2176 * t2183;
t2110 = t2226 + t2230;
t2134 = pkin(1) * t2175 + qJ(3,1);
t2170 = legFrame(1,2);
t2143 = sin(t2170);
t2146 = cos(t2170);
t2181 = cos(qJ(2,1));
t2167 = t2181 ^ 2;
t2245 = t2182 * t2184;
t2249 = t2175 * t2184;
t2290 = t2143 * qJ(3,1);
t2068 = (t2146 * t2245 - t2290) * t2167 + (t2110 * t2146 + t2143 * t2249) * t2181 + t2143 * t2134;
t2138 = t2171 * qJ(3,3);
t2238 = t2184 * t2177;
t2193 = t2138 + pkin(1) + t2238;
t2114 = 0.1e1 / t2193;
t2139 = t2173 * qJ(3,2);
t2237 = t2184 * t2179;
t2192 = t2139 + pkin(1) + t2237;
t2115 = 0.1e1 / t2192;
t2140 = t2175 * qJ(3,1);
t2236 = t2184 * t2181;
t2191 = t2140 + pkin(1) + t2236;
t2116 = 0.1e1 / t2191;
t2122 = t2146 * g(1) - t2143 * g(2);
t2200 = g(3) * t2182 + t2122 * t2176;
t2272 = t2200 * t2176;
t2214 = t2175 * t2272;
t2121 = t2145 * g(1) - t2142 * g(2);
t2201 = g(3) * t2180 + t2121 * t2174;
t2273 = t2201 * t2174;
t2217 = t2173 * t2273;
t2120 = t2144 * g(1) - t2141 * g(2);
t2202 = g(3) * t2178 + t2120 * t2172;
t2274 = t2202 * t2172;
t2220 = t2171 * t2274;
t2119 = t2143 * g(1) + t2146 * g(2);
t2303 = -g(3) * t2176 + t2122 * t2182;
t2084 = t2119 * t2175 + t2181 * t2303;
t2187 = 0.1e1 / qJ(3,1);
t2306 = t2084 * t2187;
t2118 = t2142 * g(1) + t2145 * g(2);
t2304 = -g(3) * t2174 + t2121 * t2180;
t2080 = t2118 * t2173 + t2179 * t2304;
t2186 = 0.1e1 / qJ(3,2);
t2307 = t2080 * t2186;
t2117 = t2141 * g(1) + t2144 * g(2);
t2305 = -g(3) * t2172 + t2120 * t2178;
t2076 = t2117 * t2171 + t2177 * t2305;
t2185 = 0.1e1 / qJ(3,3);
t2308 = t2076 * t2185;
t2317 = (t2064 * t2308 + t2144 * t2220) * t2114 + (t2066 * t2307 + t2145 * t2217) * t2115 + (t2068 * t2306 + t2146 * t2214) * t2116;
t2289 = t2144 * qJ(3,3);
t2063 = (-t2141 * t2247 - t2289) * t2165 + (-t2106 * t2141 + t2144 * t2255) * t2177 + t2144 * t2132;
t2288 = t2145 * qJ(3,2);
t2065 = (-t2142 * t2246 - t2288) * t2166 + (-t2108 * t2142 + t2145 * t2252) * t2179 + t2145 * t2133;
t2287 = t2146 * qJ(3,1);
t2067 = (-t2143 * t2245 - t2287) * t2167 + (-t2110 * t2143 + t2146 * t2249) * t2181 + t2146 * t2134;
t2316 = (t2063 * t2308 - t2141 * t2220) * t2114 + (t2065 * t2307 - t2142 * t2217) * t2115 + (t2067 * t2306 - t2143 * t2214) * t2116;
t2242 = t2183 * t2182;
t2089 = t2191 * t2176 - t2242;
t2275 = t2089 * t2187;
t2243 = t2183 * t2180;
t2088 = t2192 * t2174 - t2243;
t2276 = t2088 * t2186;
t2244 = t2183 * t2178;
t2087 = t2193 * t2172 - t2244;
t2277 = t2087 * t2185;
t2309 = (t2076 * t2177 * t2277 - t2202 * t2256) * t2114 + (t2080 * t2179 * t2276 - t2201 * t2253) * t2115 + (t2084 * t2181 * t2275 - t2200 * t2250) * t2116;
t2299 = -0.2e1 * t2184;
t2075 = -t2117 * t2177 + t2171 * t2305;
t2286 = t2075 * t2185;
t2079 = -t2118 * t2179 + t2173 * t2304;
t2283 = t2079 * t2186;
t2083 = -t2119 * t2181 + t2175 * t2303;
t2280 = t2083 * t2187;
t2271 = (0.2e1 * t2224 + t2232) * t2184;
t2270 = (0.2e1 * t2225 + t2231) * t2184;
t2269 = (0.2e1 * t2226 + t2230) * t2184;
t2268 = t2114 * t2172;
t2267 = t2114 * t2178;
t2266 = t2115 * t2174;
t2265 = t2115 * t2180;
t2264 = t2116 * t2176;
t2263 = t2116 * t2182;
t2259 = (qJ(3,3) + t2184) * (-qJ(3,3) + t2184);
t2258 = (qJ(3,2) + t2184) * (-qJ(3,2) + t2184);
t2257 = (qJ(3,1) + t2184) * (-qJ(3,1) + t2184);
t2235 = t2177 * pkin(2) + t2138;
t2123 = pkin(1) + t2235;
t2072 = t2120 * (-t2178 * pkin(5) + t2123 * t2172) + g(3) * (t2172 * pkin(5) + t2123 * t2178);
t2254 = t2172 * t2072;
t2234 = t2179 * pkin(2) + t2139;
t2124 = pkin(1) + t2234;
t2073 = t2121 * (-t2180 * pkin(5) + t2124 * t2174) + g(3) * (t2174 * pkin(5) + t2124 * t2180);
t2251 = t2174 * t2073;
t2233 = t2181 * pkin(2) + t2140;
t2125 = pkin(1) + t2233;
t2074 = t2122 * (-t2182 * pkin(5) + t2125 * t2176) + g(3) * (t2176 * pkin(5) + t2125 * t2182);
t2248 = t2176 * t2074;
t2241 = t2184 * t2132;
t2240 = t2184 * t2133;
t2239 = t2184 * t2134;
t2229 = qJ(3,1) * t2299;
t2228 = qJ(3,2) * t2299;
t2227 = qJ(3,3) * t2299;
t2218 = t2177 * t2274;
t2215 = t2179 * t2273;
t2212 = t2181 * t2272;
t2211 = t2141 * t2268;
t2210 = t2144 * t2268;
t2209 = t2142 * t2266;
t2208 = t2145 * t2266;
t2207 = t2143 * t2264;
t2206 = t2146 * t2264;
t2205 = t2178 * t2259;
t2204 = t2180 * t2258;
t2203 = t2182 * t2257;
t2190 = t2263 * t2303 + t2265 * t2304 + t2267 * t2305;
t2189 = t2207 * t2303 + t2209 * t2304 + t2211 * t2305;
t2188 = t2206 * t2303 + t2208 * t2304 + t2210 * t2305;
t2113 = pkin(1) * qJ(3,1) - t2175 * t2257;
t2112 = pkin(1) * qJ(3,2) - t2173 * t2258;
t2111 = pkin(1) * qJ(3,3) - t2171 * t2259;
t2104 = qJ(3,1) * t2182 + t2230 * t2175;
t2103 = qJ(3,2) * t2180 + t2231 * t2173;
t2102 = qJ(3,3) * t2178 + t2232 * t2171;
t2071 = -t2119 * t2233 + t2303 * (t2175 * pkin(2) - t2181 * qJ(3,1));
t2070 = -t2118 * t2234 + t2304 * (t2173 * pkin(2) - t2179 * qJ(3,2));
t2069 = -t2117 * t2235 + t2305 * (t2171 * pkin(2) - t2177 * qJ(3,3));
t2062 = (-t2083 * t2275 - t2182 * t2200) * t2181 * t2116 + (-t2079 * t2276 - t2180 * t2201) * t2179 * t2115 + (-t2075 * t2277 - t2178 * t2202) * t2177 * t2114;
t2061 = (t2068 * t2280 - t2146 * t2212) * t2116 + (t2066 * t2283 - t2145 * t2215) * t2115 + (t2064 * t2286 - t2144 * t2218) * t2114;
t2060 = (t2067 * t2280 + t2143 * t2212) * t2116 + (t2065 * t2283 + t2142 * t2215) * t2115 + (t2063 * t2286 + t2141 * t2218) * t2114;
t1 = [0, -t2200 * t2206 - t2201 * t2208 - t2202 * t2210, -t2188, 0, 0, 0, 0, 0, t2061, t2317, t2061, t2188, -t2317, (-t2146 * t2248 + (t2068 * t2071 - ((t2143 * t2229 + t2146 * t2203) * t2167 + (-t2113 * t2143 + t2146 * t2269) * t2181 + t2104 * t2287 + t2143 * t2239) * t2083) * t2187) * t2116 + (-t2145 * t2251 + (t2066 * t2070 - ((t2142 * t2228 + t2145 * t2204) * t2166 + (-t2112 * t2142 + t2145 * t2270) * t2179 + t2103 * t2288 + t2142 * t2240) * t2079) * t2186) * t2115 + (-t2144 * t2254 + (t2064 * t2069 - ((t2141 * t2227 + t2144 * t2205) * t2165 + (-t2111 * t2141 + t2144 * t2271) * t2177 + t2102 * t2289 + t2141 * t2241) * t2075) * t2185) * t2114, -g(1); 0, t2200 * t2207 + t2201 * t2209 + t2202 * t2211, t2189, 0, 0, 0, 0, 0, t2060, t2316, t2060, -t2189, -t2316, (t2143 * t2248 + (t2067 * t2071 - ((-t2143 * t2203 + t2146 * t2229) * t2167 + (-t2113 * t2146 - t2143 * t2269) * t2181 - t2104 * t2290 + t2146 * t2239) * t2083) * t2187) * t2116 + (t2142 * t2251 + (t2065 * t2070 - ((-t2142 * t2204 + t2145 * t2228) * t2166 + (-t2112 * t2145 - t2142 * t2270) * t2179 - t2103 * t2291 + t2145 * t2240) * t2079) * t2186) * t2115 + (t2141 * t2254 + (t2063 * t2069 - ((-t2141 * t2205 + t2144 * t2227) * t2165 + (-t2111 * t2144 - t2141 * t2271) * t2177 - t2102 * t2292 + t2144 * t2241) * t2075) * t2185) * t2114, -g(2); 0, -t2200 * t2263 - t2201 * t2265 - t2202 * t2267, -t2190, 0, 0, 0, 0, 0, t2062, -t2309, t2062, t2190, t2309, (-t2182 * t2074 + (-t2181 * t2089 * t2071 - (-t2176 * t2167 * t2257 - ((0.2e1 * t2140 + pkin(1)) * t2176 - t2242) * t2236 - qJ(3,1) * (t2134 * t2176 - t2175 * t2242)) * t2083) * t2187) * t2116 + (-t2180 * t2073 + (-t2179 * t2088 * t2070 - (-t2174 * t2166 * t2258 - ((0.2e1 * t2139 + pkin(1)) * t2174 - t2243) * t2237 - qJ(3,2) * (t2133 * t2174 - t2173 * t2243)) * t2079) * t2186) * t2115 + (-t2178 * t2072 + (-t2177 * t2087 * t2069 - (-t2172 * t2165 * t2259 - ((0.2e1 * t2138 + pkin(1)) * t2172 - t2244) * t2238 - qJ(3,3) * (t2132 * t2172 - t2171 * t2244)) * t2075) * t2185) * t2114, -g(3);];
tau_reg  = t1;
