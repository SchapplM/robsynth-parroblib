% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR2G3A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR2G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:51
% EndTime: 2020-03-09 21:11:57
% DurationCPUTime: 5.62s
% Computational Cost: add. (11547->466), mult. (27675->852), div. (10512->11), fcn. (30786->42), ass. (0->346)
t224 = xDP(2);
t217 = cos(qJ(3,1));
t190 = 0.1e1 / t217 ^ 2;
t209 = sin(qJ(2,1));
t181 = 0.1e1 / t209;
t228 = 0.1e1 / pkin(2);
t230 = 0.1e1 / pkin(1);
t295 = t228 * t230;
t247 = t181 * t295;
t241 = t190 * t247;
t198 = legFrame(1,2);
t174 = sin(t198);
t177 = cos(t198);
t210 = sin(qJ(1,1));
t218 = cos(qJ(2,1));
t219 = cos(qJ(1,1));
t132 = t209 * t210 - t218 * t219;
t188 = t217 ^ 2;
t286 = pkin(2) * t132 * t188;
t208 = sin(qJ(3,1));
t374 = pkin(1) * t218;
t289 = t208 * t374;
t321 = t177 * t208;
t373 = pkin(1) * t219;
t81 = -t174 * t286 + (-pkin(2) * t321 + t174 * t373) * t217 - t177 * t289;
t75 = t81 * t224 * t241;
t225 = xDP(1);
t324 = t174 * t208;
t84 = t177 * t286 + (-pkin(2) * t324 - t177 * t373) * t217 - t174 * t289;
t78 = t84 * t225 * t241;
t117 = pkin(2) * (t209 * t219 + t210 * t218) * t217 + t210 * pkin(1);
t223 = xDP(3);
t189 = 0.1e1 / t217;
t242 = t189 * t247;
t99 = t117 * t223 * t242;
t51 = t78 + t75 + t99;
t330 = t132 * t217;
t104 = t174 * t330 + t321;
t105 = -t177 * t330 + t324;
t193 = qJ(1,1) + qJ(2,1);
t161 = sin(t193);
t318 = t181 * t230;
t248 = t189 * t318;
t296 = t223 * t230;
t69 = -t161 * t181 * t296 + (t104 * t224 + t105 * t225) * t248;
t369 = -t51 - t69;
t33 = t369 ^ 2;
t214 = cos(qJ(3,2));
t187 = 0.1e1 / t214 ^ 2;
t206 = sin(qJ(2,2));
t180 = 0.1e1 / t206;
t249 = t180 * t295;
t243 = t187 * t249;
t197 = legFrame(2,2);
t173 = sin(t197);
t176 = cos(t197);
t207 = sin(qJ(1,2));
t215 = cos(qJ(2,2));
t216 = cos(qJ(1,2));
t131 = t206 * t207 - t215 * t216;
t185 = t214 ^ 2;
t287 = pkin(2) * t131 * t185;
t205 = sin(qJ(3,2));
t376 = pkin(1) * t215;
t290 = t205 * t376;
t322 = t176 * t205;
t375 = pkin(1) * t216;
t80 = -t173 * t287 + (-pkin(2) * t322 + t173 * t375) * t214 - t176 * t290;
t74 = t80 * t224 * t243;
t325 = t173 * t205;
t83 = t176 * t287 + (-pkin(2) * t325 - t176 * t375) * t214 - t173 * t290;
t77 = t83 * t225 * t243;
t116 = pkin(2) * (t206 * t216 + t207 * t215) * t214 + t207 * pkin(1);
t186 = 0.1e1 / t214;
t244 = t186 * t249;
t98 = t116 * t223 * t244;
t50 = t77 + t74 + t98;
t331 = t131 * t214;
t102 = t173 * t331 + t322;
t103 = -t176 * t331 + t325;
t192 = qJ(1,2) + qJ(2,2);
t160 = sin(t192);
t319 = t180 * t230;
t250 = t186 * t319;
t68 = -t160 * t180 * t296 + (t102 * t224 + t103 * t225) * t250;
t370 = -t50 - t68;
t32 = t370 ^ 2;
t211 = cos(qJ(3,3));
t184 = 0.1e1 / t211 ^ 2;
t203 = sin(qJ(2,3));
t179 = 0.1e1 / t203;
t251 = t179 * t295;
t245 = t184 * t251;
t196 = legFrame(3,2);
t172 = sin(t196);
t175 = cos(t196);
t204 = sin(qJ(1,3));
t212 = cos(qJ(2,3));
t213 = cos(qJ(1,3));
t130 = t203 * t204 - t212 * t213;
t182 = t211 ^ 2;
t288 = pkin(2) * t130 * t182;
t202 = sin(qJ(3,3));
t378 = pkin(1) * t212;
t291 = t202 * t378;
t323 = t175 * t202;
t377 = pkin(1) * t213;
t79 = -t172 * t288 + (-pkin(2) * t323 + t172 * t377) * t211 - t175 * t291;
t73 = t79 * t224 * t245;
t326 = t172 * t202;
t82 = t175 * t288 + (-pkin(2) * t326 - t175 * t377) * t211 - t172 * t291;
t76 = t82 * t225 * t245;
t115 = pkin(2) * (t203 * t213 + t204 * t212) * t211 + t204 * pkin(1);
t183 = 0.1e1 / t211;
t246 = t183 * t251;
t97 = t115 * t223 * t246;
t49 = t76 + t73 + t97;
t332 = t130 * t211;
t100 = t172 * t332 + t323;
t101 = -t175 * t332 + t326;
t191 = qJ(1,3) + qJ(2,3);
t159 = sin(t191);
t320 = t179 * t230;
t252 = t183 * t320;
t67 = -t159 * t179 * t296 + (t100 * t224 + t101 * t225) * t252;
t371 = -t49 - t67;
t31 = t371 ^ 2;
t135 = g(1) * t177 - g(2) * t174;
t195 = mrSges(3,3) - mrSges(2,2);
t222 = mrSges(2,1) * g(3);
t393 = -t135 * t195 + t222;
t134 = g(1) * t176 - g(2) * t173;
t392 = -t134 * t195 + t222;
t133 = g(1) * t175 - g(2) * t172;
t391 = -t133 * t195 + t222;
t390 = 0.2e1 * pkin(1);
t64 = t67 ^ 2;
t389 = t64 * t378;
t65 = t68 ^ 2;
t388 = t65 * t376;
t66 = t69 ^ 2;
t387 = t66 * t374;
t171 = g(3) * t195;
t386 = mrSges(2,1) * t133 + t171;
t385 = mrSges(2,1) * t134 + t171;
t384 = mrSges(2,1) * t135 + t171;
t383 = g(3) * mrSges(1,2);
t381 = pkin(1) * t203;
t380 = pkin(1) * t206;
t379 = pkin(1) * t209;
t372 = Ifges(3,1) + Ifges(2,3);
t368 = mrSges(3,2) * t133;
t367 = mrSges(3,2) * t134;
t366 = mrSges(3,2) * t135;
t365 = mrSges(3,2) * t202;
t364 = mrSges(3,2) * t205;
t363 = mrSges(3,2) * t208;
t362 = Ifges(3,4) * t202;
t361 = Ifges(3,4) * t205;
t360 = Ifges(3,4) * t208;
t317 = t183 * t228;
t112 = (t172 * t225 + t175 * t224) * t317;
t359 = t112 * t371;
t315 = t186 * t228;
t113 = (t173 * t225 + t176 * t224) * t315;
t358 = t113 * t370;
t313 = t189 * t228;
t114 = (t174 * t225 + t177 * t224) * t313;
t357 = t114 * t369;
t109 = t112 ^ 2;
t356 = (t109 / 0.2e1 + (t67 + t49 / 0.2e1) * t49) * t203;
t110 = t113 ^ 2;
t355 = (t110 / 0.2e1 + (t68 + t50 / 0.2e1) * t50) * t206;
t111 = t114 ^ 2;
t354 = (t111 / 0.2e1 + (t69 + t51 / 0.2e1) * t51) * t209;
t353 = t182 * t371;
t118 = t211 * (-mrSges(3,2) * t381 + Ifges(3,6)) - t202 * (mrSges(3,1) * t381 - Ifges(3,5));
t142 = Ifges(3,5) * t202 + Ifges(3,6) * t211;
t267 = t115 * t317;
t352 = (-t118 * t159 + t142 * t267) * t252;
t351 = t185 * t370;
t119 = t214 * (-mrSges(3,2) * t380 + Ifges(3,6)) - t205 * (mrSges(3,1) * t380 - Ifges(3,5));
t143 = Ifges(3,5) * t205 + Ifges(3,6) * t214;
t266 = t116 * t315;
t350 = (-t119 * t160 + t143 * t266) * t250;
t349 = t188 * t369;
t120 = t217 * (-mrSges(3,2) * t379 + Ifges(3,6)) - t208 * (mrSges(3,1) * t379 - Ifges(3,5));
t144 = Ifges(3,5) * t208 + Ifges(3,6) * t217;
t265 = t117 * t313;
t348 = (-t120 * t161 + t144 * t265) * t248;
t347 = t100 * t183;
t346 = t101 * t183;
t345 = t102 * t186;
t344 = t103 * t186;
t343 = t104 * t189;
t342 = t105 * t189;
t341 = t109 * t202;
t340 = t110 * t205;
t339 = t111 * t208;
t338 = t112 * t212;
t337 = t113 * t215;
t336 = t114 * t218;
t335 = t118 * t183;
t334 = t119 * t186;
t333 = t120 * t189;
t316 = t184 * t228;
t314 = t187 * t228;
t312 = t190 * t228;
t194 = Ifges(3,2) - Ifges(3,1);
t311 = t194 * t182;
t310 = t194 * t185;
t309 = t194 * t188;
t308 = t194 * t202;
t307 = t194 * t205;
t306 = t194 * t208;
t305 = t195 * t212;
t304 = t195 * t215;
t303 = t195 * t218;
t302 = t202 * t203;
t301 = t205 * t206;
t300 = t208 * t209;
t299 = t211 * t212;
t298 = t214 * t215;
t297 = t217 * t218;
t294 = 0.2e1 * t362;
t293 = 0.2e1 * t361;
t292 = 0.2e1 * t360;
t282 = t371 * t338;
t281 = t370 * t337;
t280 = t369 * t336;
t279 = t79 * t316;
t278 = t82 * t316;
t136 = (-mrSges(2,1) + t365) * t378;
t149 = t195 * t381;
t156 = mrSges(3,1) * t378;
t240 = t311 + t372;
t88 = (t156 + t294) * t211 - t136 + t149 + t240;
t277 = t88 * t316;
t276 = t80 * t314;
t275 = t83 * t314;
t137 = (-mrSges(2,1) + t364) * t376;
t150 = t195 * t380;
t157 = mrSges(3,1) * t376;
t239 = t310 + t372;
t89 = (t157 + t293) * t214 - t137 + t150 + t239;
t274 = t89 * t314;
t273 = t81 * t312;
t272 = t84 * t312;
t138 = (-mrSges(2,1) + t363) * t374;
t151 = t195 * t379;
t158 = mrSges(3,1) * t374;
t238 = t309 + t372;
t90 = (t158 + t292) * t217 - t138 + t151 + t238;
t271 = t90 * t312;
t270 = t112 * t302;
t269 = t113 * t301;
t268 = t114 * t300;
t121 = t211 * t294 + t240;
t264 = t121 * t316;
t122 = t214 * t293 + t239;
t263 = t122 * t314;
t123 = t217 * t292 + t238;
t262 = t123 * t312;
t261 = t142 * t316;
t260 = t143 * t314;
t259 = t144 * t312;
t258 = t172 * t317;
t257 = t173 * t315;
t256 = t174 * t313;
t255 = t175 * t317;
t254 = t176 * t315;
t253 = t177 * t313;
t237 = mrSges(3,1) * t211 - t365;
t236 = mrSges(3,1) * t214 - t364;
t235 = mrSges(3,1) * t217 - t363;
t226 = m(2) + m(3);
t234 = pkin(1) ^ 2 * t226 + Ifges(1,3) + t372;
t233 = (Ifges(3,5) * t109 + 0.2e1 * t308 * t359) * t211 - (0.4e1 * t182 - 0.2e1) * Ifges(3,4) * t359;
t232 = (Ifges(3,5) * t110 + 0.2e1 * t307 * t358) * t214 - (0.4e1 * t185 - 0.2e1) * Ifges(3,4) * t358;
t231 = (Ifges(3,5) * t111 + 0.2e1 * t306 * t357) * t217 - (0.4e1 * t188 - 0.2e1) * Ifges(3,4) * t357;
t227 = pkin(2) ^ 2;
t221 = mrSges(3,1) * g(3);
t220 = mrSges(3,2) * g(3);
t201 = xDDP(1);
t200 = xDDP(2);
t199 = xDDP(3);
t170 = -qJ(3,1) + t193;
t169 = qJ(3,1) + t193;
t168 = -qJ(3,2) + t192;
t167 = qJ(3,2) + t192;
t166 = -qJ(3,3) + t191;
t165 = qJ(3,3) + t191;
t164 = cos(t193);
t163 = cos(t192);
t162 = cos(t191);
t155 = pkin(1) * t226 + mrSges(1,1);
t148 = t155 * g(3);
t126 = mrSges(3,1) * t135;
t125 = mrSges(3,1) * t134;
t124 = mrSges(3,1) * t133;
t87 = t309 + 0.2e1 * (t158 + t360) * t217 - 0.2e1 * t138 + 0.2e1 * t151 + t234;
t86 = t310 + 0.2e1 * (t157 + t361) * t214 - 0.2e1 * t137 + 0.2e1 * t150 + t234;
t85 = t311 + 0.2e1 * (t156 + t362) * t211 - 0.2e1 * t136 + 0.2e1 * t149 + t234;
t63 = (t123 * t265 - t161 * t90) * t318;
t62 = (t122 * t266 - t160 * t89) * t319;
t61 = (t121 * t267 - t159 * t88) * t320;
t60 = (-t161 * t87 + t265 * t90) * t318;
t59 = (-t160 * t86 + t266 * t89) * t319;
t58 = (-t159 * t85 + t267 * t88) * t320;
t57 = Ifges(3,3) * t256 + (t105 * t333 + t259 * t84) * t318;
t56 = Ifges(3,3) * t253 + (t104 * t333 + t259 * t81) * t318;
t55 = Ifges(3,3) * t257 + (t103 * t334 + t260 * t83) * t319;
t54 = Ifges(3,3) * t254 + (t102 * t334 + t260 * t80) * t319;
t53 = Ifges(3,3) * t258 + (t101 * t335 + t261 * t82) * t320;
t52 = Ifges(3,3) * t255 + (t100 * t335 + t261 * t79) * t320;
t48 = t144 * t256 + (t262 * t84 + t342 * t90) * t318;
t47 = t144 * t253 + (t262 * t81 + t343 * t90) * t318;
t46 = t143 * t257 + (t263 * t83 + t344 * t89) * t319;
t45 = t143 * t254 + (t263 * t80 + t345 * t89) * t319;
t44 = t142 * t258 + (t264 * t82 + t346 * t88) * t320;
t43 = t142 * t255 + (t264 * t79 + t347 * t88) * t320;
t42 = t120 * t256 + (t271 * t84 + t342 * t87) * t318;
t41 = t120 * t253 + (t271 * t81 + t343 * t87) * t318;
t40 = t119 * t257 + (t274 * t83 + t344 * t86) * t319;
t39 = t119 * t254 + (t274 * t80 + t345 * t86) * t319;
t38 = t118 * t258 + (t277 * t82 + t346 * t85) * t320;
t37 = t118 * t255 + (t277 * t79 + t347 * t85) * t320;
t30 = t78 / 0.2e1 + t75 / 0.2e1 + t99 / 0.2e1 + t69;
t29 = t77 / 0.2e1 + t74 / 0.2e1 + t98 / 0.2e1 + t68;
t28 = t76 / 0.2e1 + t73 / 0.2e1 + t97 / 0.2e1 + t67;
t15 = (-t387 + (-t111 * t189 - t217 * t33) * pkin(2)) * t318;
t14 = (-t388 + (-t110 * t186 - t214 * t32) * pkin(2)) * t319;
t13 = (-t389 + (-t109 * t183 - t211 * t31) * pkin(2)) * t320;
t12 = (-t227 * t349 + (pkin(1) * t69 + (0.2e1 * t297 * t30 - t268) * pkin(2)) * pkin(1)) * t69 * t242 + (-pkin(2) * t349 + (-t297 * t369 - t268) * pkin(1)) * t51 * t248 + ((pkin(1) * t300 * t369 + pkin(2) * t114) * t217 + pkin(1) * t336) * t190 * t114 * t318;
t11 = (-t227 * t351 + (pkin(1) * t68 + (0.2e1 * t29 * t298 - t269) * pkin(2)) * pkin(1)) * t68 * t244 + (-pkin(2) * t351 + (-t298 * t370 - t269) * pkin(1)) * t50 * t250 + ((pkin(1) * t301 * t370 + pkin(2) * t113) * t214 + pkin(1) * t337) * t187 * t113 * t319;
t10 = (-t227 * t353 + (pkin(1) * t67 + (0.2e1 * t28 * t299 - t270) * pkin(2)) * pkin(1)) * t67 * t246 + (-pkin(2) * t353 + (-t299 * t371 - t270) * pkin(1)) * t49 * t252 + ((pkin(1) * t302 * t371 + pkin(2) * t112) * t211 + pkin(1) * t338) * t184 * t112 * t320;
t9 = -t120 * t15 - t144 * t12 + Ifges(3,3) * t189 * t339 + (mrSges(3,2) * t387 + t33 * t306) * t217 + mrSges(3,1) * t66 * t289 - (g(1) * t174 + g(2) * t177) * t235 + (-g(3) * t161 + t135 * t164) * (mrSges(3,1) * t208 + mrSges(3,2) * t217) + (-0.2e1 * t188 + 0.1e1) * t33 * Ifges(3,4);
t8 = -t119 * t14 - t143 * t11 + Ifges(3,3) * t186 * t340 + (mrSges(3,2) * t388 + t32 * t307) * t214 + mrSges(3,1) * t65 * t290 - (g(1) * t173 + g(2) * t176) * t236 + (-g(3) * t160 + t134 * t163) * (mrSges(3,1) * t205 + mrSges(3,2) * t214) + (-0.2e1 * t185 + 0.1e1) * t32 * Ifges(3,4);
t7 = -t118 * t13 - t142 * t10 + Ifges(3,3) * t183 * t341 + (mrSges(3,2) * t389 + t31 * t308) * t211 + mrSges(3,1) * t64 * t291 - (g(1) * t172 + g(2) * t175) * t237 + (-g(3) * t159 + t133 * t162) * (mrSges(3,1) * t202 + mrSges(3,2) * t211) + (-0.2e1 * t182 + 0.1e1) * t31 * Ifges(3,4);
t6 = -t90 * t15 - t123 * t12 + (g(3) * t235 + t393) * t164 + t161 * (t135 * t235 + t384) + (t144 * t189 - Ifges(3,6)) * t339 + (-t303 + (mrSges(2,1) + t235) * t209) * t66 * pkin(1) + t231;
t5 = -t89 * t14 - t122 * t11 + (g(3) * t236 + t392) * t163 + t160 * (t134 * t236 + t385) + (t143 * t186 - Ifges(3,6)) * t340 + (-t304 + (mrSges(2,1) + t236) * t206) * t65 * pkin(1) + t232;
t4 = -t88 * t13 - t121 * t10 + (g(3) * t237 + t391) * t162 + t159 * (t133 * t237 + t386) + (t142 * t183 - Ifges(3,6)) * t341 + (-t305 + (mrSges(2,1) + t237) * t203) * t64 * pkin(1) + t233;
t3 = -t87 * t15 - t90 * t12 + (t221 - t366) * cos(t170) / 0.2e1 + (t126 + t220) * sin(t170) / 0.2e1 + (t221 + t366) * cos(t169) / 0.2e1 + (t126 - t220) * sin(t169) / 0.2e1 + (t135 * t155 - t383) * t210 + t384 * t161 + (mrSges(1,2) * t135 + t148) * t219 + (-Ifges(3,6) + t333) * t339 + ((-mrSges(3,1) * t354 + mrSges(3,2) * t280) * t217 + (mrSges(3,1) * t280 + mrSges(3,2) * t354) * t208 + (-mrSges(2,1) * t209 + t303) * t51 * t30) * t390 + t231 + t393 * t164;
t2 = -t86 * t14 - t89 * t11 + (t221 - t367) * cos(t168) / 0.2e1 + (t125 + t220) * sin(t168) / 0.2e1 + (t221 + t367) * cos(t167) / 0.2e1 + (t125 - t220) * sin(t167) / 0.2e1 + (t134 * t155 - t383) * t207 + t385 * t160 + (mrSges(1,2) * t134 + t148) * t216 + (-Ifges(3,6) + t334) * t340 + ((-mrSges(3,1) * t355 + mrSges(3,2) * t281) * t214 + (mrSges(3,1) * t281 + mrSges(3,2) * t355) * t205 + (-mrSges(2,1) * t206 + t304) * t50 * t29) * t390 + t232 + t392 * t163;
t1 = -t85 * t13 - t88 * t10 + (t221 - t368) * cos(t166) / 0.2e1 + (t124 + t220) * sin(t166) / 0.2e1 + (t221 + t368) * cos(t165) / 0.2e1 + (t124 - t220) * sin(t165) / 0.2e1 + (t133 * t155 - t383) * t204 + t386 * t159 + (mrSges(1,2) * t133 + t148) * t213 + (-Ifges(3,6) + t335) * t341 + ((-mrSges(3,1) * t356 + mrSges(3,2) * t282) * t211 + (mrSges(3,1) * t282 + mrSges(3,2) * t356) * t202 + (-mrSges(2,1) * t203 + t305) * t49 * t28) * t390 + t233 + t391 * t162;
t16 = [(-g(1) + t201) * m(4) + ((t177 * t200 * t57 + (t201 * t57 + t9) * t174) * t189 + (t176 * t200 * t55 + (t201 * t55 + t8) * t173) * t186 + (t175 * t200 * t53 + (t201 * t53 + t7) * t172) * t183) * t228 + (((t272 * t48 + t342 * t42) * t201 + (t273 * t48 + t343 * t42) * t200 + (-t161 * t42 + t265 * t48) * t199 + t3 * t342 + t6 * t272) * t181 + ((t275 * t46 + t344 * t40) * t201 + (t276 * t46 + t345 * t40) * t200 + (-t160 * t40 + t266 * t46) * t199 + t2 * t344 + t5 * t275) * t180 + ((t278 * t44 + t346 * t38) * t201 + (t279 * t44 + t347 * t38) * t200 + (-t159 * t38 + t267 * t44) * t199 + t1 * t346 + t4 * t278) * t179) * t230; (-g(2) + t200) * m(4) + ((t174 * t201 * t56 + (t200 * t56 + t9) * t177) * t189 + (t173 * t201 * t54 + (t200 * t54 + t8) * t176) * t186 + (t172 * t201 * t52 + (t200 * t52 + t7) * t175) * t183) * t228 + (((t272 * t47 + t342 * t41) * t201 + (t273 * t47 + t343 * t41) * t200 + (-t161 * t41 + t265 * t47) * t199 + t3 * t343 + t6 * t273) * t181 + ((t275 * t45 + t344 * t39) * t201 + (t276 * t45 + t345 * t39) * t200 + (-t160 * t39 + t266 * t45) * t199 + t2 * t345 + t5 * t276) * t180 + ((t278 * t43 + t346 * t37) * t201 + (t279 * t43 + t347 * t37) * t200 + (-t159 * t37 + t267 * t43) * t199 + t1 * t347 + t4 * t279) * t179) * t230; (-g(3) + t199) * m(4) + ((t172 * t352 + t173 * t350 + t174 * t348) * t201 + (t175 * t352 + t176 * t350 + t177 * t348) * t200) * t228 + (((t272 * t63 + t342 * t60) * t201 + (t273 * t63 + t343 * t60) * t200 + (-t161 * t60 + t265 * t63) * t199 - t161 * t3 + t6 * t265) * t181 + ((t275 * t62 + t344 * t59) * t201 + (t276 * t62 + t345 * t59) * t200 + (-t160 * t59 + t266 * t62) * t199 - t160 * t2 + t5 * t266) * t180 + ((t278 * t61 + t346 * t58) * t201 + (t279 * t61 + t347 * t58) * t200 + (-t159 * t58 + t267 * t61) * t199 - t159 * t1 + t4 * t267) * t179) * t230;];
tauX  = t16;
