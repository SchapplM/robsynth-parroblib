% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR2G2A0
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
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR2G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:44
% EndTime: 2020-03-09 21:08:49
% DurationCPUTime: 5.64s
% Computational Cost: add. (14445->478), mult. (27261->858), div. (10098->14), fcn. (29544->48), ass. (0->363)
t228 = xDP(2);
t220 = cos(qJ(3,1));
t193 = 0.1e1 / t220 ^ 2;
t212 = sin(qJ(2,1));
t184 = 0.1e1 / t212;
t232 = 0.1e1 / pkin(2);
t234 = 0.1e1 / pkin(1);
t306 = t232 * t234;
t252 = t184 * t306;
t248 = t193 * t252;
t201 = legFrame(1,2);
t178 = sin(t201);
t181 = cos(t201);
t213 = sin(qJ(1,1));
t221 = cos(qJ(2,1));
t222 = cos(qJ(1,1));
t126 = t212 * t222 + t213 * t221;
t191 = t220 ^ 2;
t291 = pkin(2) * t126 * t191;
t211 = sin(qJ(3,1));
t397 = pkin(1) * t221;
t294 = t211 * t397;
t333 = t181 * t211;
t400 = pkin(1) * t213;
t83 = t178 * t291 + (-pkin(2) * t333 + t178 * t400) * t220 - t181 * t294;
t77 = t83 * t228 * t248;
t229 = xDP(1);
t336 = t178 * t211;
t84 = -t181 * t291 + (-pkin(2) * t336 - t181 * t400) * t220 - t178 * t294;
t78 = t84 * t229 * t248;
t227 = xDP(3);
t251 = t227 * t306;
t304 = qJ(2,1) + qJ(3,1);
t173 = qJ(1,1) + t304;
t157 = cos(t173);
t305 = qJ(2,1) - qJ(3,1);
t174 = qJ(1,1) + t305;
t158 = cos(t174);
t407 = -0.2e1 * pkin(1);
t348 = (t222 * t407 + (-t157 - t158) * pkin(2)) / (sin(t304) + sin(t305));
t99 = t251 * t348;
t51 = t78 + t77 + t99;
t345 = t126 * t220;
t104 = -t178 * t345 + t333;
t105 = t181 * t345 + t336;
t192 = 0.1e1 / t220;
t330 = t184 * t234;
t253 = t192 * t330;
t307 = t227 * t234;
t196 = qJ(1,1) + qJ(2,1);
t168 = cos(t196);
t340 = t168 * t184;
t69 = t307 * t340 + (t104 * t228 + t105 * t229) * t253;
t393 = -t51 - t69;
t33 = t393 ^ 2;
t217 = cos(qJ(3,2));
t190 = 0.1e1 / t217 ^ 2;
t209 = sin(qJ(2,2));
t183 = 0.1e1 / t209;
t254 = t183 * t306;
t249 = t190 * t254;
t200 = legFrame(2,2);
t177 = sin(t200);
t180 = cos(t200);
t210 = sin(qJ(1,2));
t218 = cos(qJ(2,2));
t219 = cos(qJ(1,2));
t125 = t209 * t219 + t210 * t218;
t188 = t217 ^ 2;
t292 = pkin(2) * t125 * t188;
t208 = sin(qJ(3,2));
t398 = pkin(1) * t218;
t295 = t208 * t398;
t334 = t180 * t208;
t402 = pkin(1) * t210;
t81 = t177 * t292 + (-pkin(2) * t334 + t177 * t402) * t217 - t180 * t295;
t75 = t81 * t228 * t249;
t337 = t177 * t208;
t82 = -t180 * t292 + (-pkin(2) * t337 - t180 * t402) * t217 - t177 * t295;
t76 = t82 * t229 * t249;
t302 = qJ(2,2) + qJ(3,2);
t171 = qJ(1,2) + t302;
t155 = cos(t171);
t303 = qJ(2,2) - qJ(3,2);
t172 = qJ(1,2) + t303;
t156 = cos(t172);
t349 = (t219 * t407 + (-t155 - t156) * pkin(2)) / (sin(t302) + sin(t303));
t98 = t251 * t349;
t50 = t76 + t75 + t98;
t346 = t125 * t217;
t102 = -t177 * t346 + t334;
t103 = t180 * t346 + t337;
t189 = 0.1e1 / t217;
t331 = t183 * t234;
t255 = t189 * t331;
t195 = qJ(1,2) + qJ(2,2);
t167 = cos(t195);
t342 = t167 * t183;
t68 = t307 * t342 + (t102 * t228 + t103 * t229) * t255;
t394 = -t50 - t68;
t32 = t394 ^ 2;
t214 = cos(qJ(3,3));
t187 = 0.1e1 / t214 ^ 2;
t206 = sin(qJ(2,3));
t182 = 0.1e1 / t206;
t256 = t182 * t306;
t250 = t187 * t256;
t199 = legFrame(3,2);
t176 = sin(t199);
t179 = cos(t199);
t207 = sin(qJ(1,3));
t215 = cos(qJ(2,3));
t216 = cos(qJ(1,3));
t124 = t206 * t216 + t207 * t215;
t185 = t214 ^ 2;
t293 = pkin(2) * t124 * t185;
t205 = sin(qJ(3,3));
t399 = pkin(1) * t215;
t296 = t205 * t399;
t335 = t179 * t205;
t404 = pkin(1) * t207;
t79 = t176 * t293 + (-pkin(2) * t335 + t176 * t404) * t214 - t179 * t296;
t73 = t79 * t228 * t250;
t338 = t176 * t205;
t80 = -t179 * t293 + (-pkin(2) * t338 - t179 * t404) * t214 - t176 * t296;
t74 = t80 * t229 * t250;
t300 = qJ(2,3) + qJ(3,3);
t169 = qJ(1,3) + t300;
t153 = cos(t169);
t301 = qJ(2,3) - qJ(3,3);
t170 = qJ(1,3) + t301;
t154 = cos(t170);
t350 = (t216 * t407 + (-t153 - t154) * pkin(2)) / (sin(t300) + sin(t301));
t97 = t251 * t350;
t49 = t74 + t73 + t97;
t347 = t124 * t214;
t100 = -t176 * t347 + t335;
t101 = t179 * t347 + t338;
t186 = 0.1e1 / t214;
t332 = t182 * t234;
t257 = t186 * t332;
t194 = qJ(1,3) + qJ(2,3);
t166 = cos(t194);
t344 = t166 * t182;
t67 = t307 * t344 + (t100 * t228 + t101 * t229) * t257;
t395 = -t49 - t67;
t31 = t395 ^ 2;
t411 = 0.2e1 * pkin(1);
t64 = t67 ^ 2;
t410 = t64 * t399;
t65 = t68 ^ 2;
t409 = t65 * t398;
t66 = t69 ^ 2;
t408 = t66 * t397;
t381 = mrSges(3,2) * t211;
t242 = mrSges(3,1) * t220 - t381;
t382 = mrSges(3,2) * t208;
t243 = mrSges(3,1) * t217 - t382;
t383 = mrSges(3,2) * t205;
t244 = mrSges(3,1) * t214 - t383;
t406 = -2 * mrSges(2,1);
t405 = pkin(1) * t206;
t403 = pkin(1) * t209;
t401 = pkin(1) * t212;
t198 = mrSges(3,3) - mrSges(2,2);
t175 = (t198 * g(3));
t396 = Ifges(3,1) + Ifges(2,3);
t127 = g(1) * t179 - g(2) * t176;
t392 = mrSges(3,1) * t127;
t128 = g(1) * t180 - g(2) * t177;
t391 = mrSges(3,1) * t128;
t129 = g(1) * t181 - g(2) * t178;
t390 = mrSges(3,1) * t129;
t386 = mrSges(3,2) * t127;
t385 = mrSges(3,2) * t128;
t384 = mrSges(3,2) * t129;
t380 = Ifges(3,4) * t205;
t379 = Ifges(3,4) * t208;
t378 = Ifges(3,4) * t211;
t329 = t186 * t232;
t112 = (t176 * t229 + t179 * t228) * t329;
t377 = t112 * t395;
t327 = t189 * t232;
t113 = (t177 * t229 + t180 * t228) * t327;
t376 = t113 * t394;
t325 = t192 * t232;
t114 = (t178 * t229 + t181 * t228) * t325;
t375 = t114 * t393;
t109 = t112 ^ 2;
t374 = (t109 / 0.2e1 + (t67 + t49 / 0.2e1) * t49) * t206;
t110 = t113 ^ 2;
t373 = (t110 / 0.2e1 + (t68 + t50 / 0.2e1) * t50) * t209;
t111 = t114 ^ 2;
t372 = (t111 / 0.2e1 + (t69 + t51 / 0.2e1) * t51) * t212;
t371 = t185 * t395;
t115 = (-mrSges(3,2) * t405 + Ifges(3,6)) * t214 - t205 * (mrSges(3,1) * t405 - Ifges(3,5));
t142 = Ifges(3,5) * t205 + Ifges(3,6) * t214;
t272 = t232 * t350;
t370 = t186 * (t115 * t344 + t142 * t272) * t234;
t369 = t188 * t394;
t116 = (-mrSges(3,2) * t403 + Ifges(3,6)) * t217 - t208 * (mrSges(3,1) * t403 - Ifges(3,5));
t143 = Ifges(3,5) * t208 + Ifges(3,6) * t217;
t271 = t232 * t349;
t368 = t189 * (t116 * t342 + t143 * t271) * t234;
t367 = t191 * t393;
t117 = (-mrSges(3,2) * t401 + Ifges(3,6)) * t220 - t211 * (mrSges(3,1) * t401 - Ifges(3,5));
t144 = Ifges(3,5) * t211 + Ifges(3,6) * t220;
t270 = t232 * t348;
t366 = t192 * (t117 * t340 + t144 * t270) * t234;
t365 = t100 * t186;
t364 = t101 * t186;
t363 = t102 * t189;
t362 = t103 * t189;
t361 = t104 * t192;
t360 = t105 * t192;
t359 = t109 * t205;
t358 = t110 * t208;
t357 = t111 * t211;
t356 = t112 * t215;
t355 = t113 * t218;
t354 = t114 * t221;
t353 = t115 * t186;
t352 = t116 * t189;
t351 = t117 * t192;
t202 = xDDP(3);
t343 = t166 * t202;
t341 = t167 * t202;
t339 = t168 * t202;
t328 = t187 * t232;
t326 = t190 * t232;
t324 = t193 * t232;
t197 = Ifges(3,2) - Ifges(3,1);
t323 = t197 * t185;
t322 = t197 * t188;
t321 = t197 * t191;
t320 = t197 * t205;
t319 = t197 * t208;
t318 = t197 * t211;
t317 = t198 * t215;
t316 = t198 * t218;
t315 = t198 * t221;
t314 = t202 * t232;
t313 = t205 * t206;
t312 = t208 * t209;
t311 = t211 * t212;
t310 = t214 * t215;
t309 = t217 * t218;
t308 = t220 * t221;
t299 = 0.2e1 * t380;
t298 = 0.2e1 * t379;
t297 = 0.2e1 * t378;
t287 = t395 * t356;
t286 = t394 * t355;
t285 = t393 * t354;
t284 = t79 * t328;
t283 = t80 * t328;
t136 = (-mrSges(2,1) + t383) * t399;
t149 = t198 * t405;
t160 = mrSges(3,1) * t399;
t247 = t323 + t396;
t88 = (t160 + t299) * t214 - t136 + t149 + t247;
t282 = t88 * t328;
t281 = t81 * t326;
t280 = t82 * t326;
t137 = (-mrSges(2,1) + t382) * t398;
t150 = t198 * t403;
t161 = mrSges(3,1) * t398;
t246 = t322 + t396;
t89 = (t161 + t298) * t217 - t137 + t150 + t246;
t279 = t89 * t326;
t278 = t83 * t324;
t277 = t84 * t324;
t138 = (-mrSges(2,1) + t381) * t397;
t151 = t198 * t401;
t162 = mrSges(3,1) * t397;
t245 = t321 + t396;
t90 = (t162 + t297) * t220 - t138 + t151 + t245;
t276 = t90 * t324;
t275 = t112 * t313;
t274 = t113 * t312;
t273 = t114 * t311;
t121 = t214 * t299 + t247;
t269 = t121 * t328;
t122 = t217 * t298 + t246;
t268 = t122 * t326;
t123 = t220 * t297 + t245;
t267 = t123 * t324;
t266 = t142 * t328;
t265 = t143 * t326;
t264 = t144 * t324;
t263 = t176 * t329;
t262 = t177 * t327;
t261 = t178 * t325;
t260 = t179 * t329;
t259 = t180 * t327;
t258 = t181 * t325;
t230 = m(2) + m(3);
t241 = pkin(1) ^ 2 * t230 + Ifges(1,3) + t396;
t240 = (t109 * Ifges(3,5) + 0.2e1 * t320 * t377) * t214 - (0.4e1 * t185 - 0.2e1) * Ifges(3,4) * t377;
t239 = (t110 * Ifges(3,5) + 0.2e1 * t319 * t376) * t217 - (0.4e1 * t188 - 0.2e1) * Ifges(3,4) * t376;
t238 = (t111 * Ifges(3,5) + 0.2e1 * t318 * t375) * t220 - (0.4e1 * t191 - 0.2e1) * Ifges(3,4) * t375;
t237 = mrSges(2,1) + t244;
t236 = mrSges(2,1) + t243;
t235 = mrSges(2,1) + t242;
t231 = pkin(2) ^ 2;
t226 = mrSges(2,1) * g(3);
t225 = mrSges(3,1) * g(3);
t224 = mrSges(1,2) * g(3);
t223 = mrSges(3,2) * g(3);
t204 = xDDP(1);
t203 = xDDP(2);
t165 = sin(t196);
t164 = sin(t195);
t163 = sin(t194);
t159 = pkin(1) * t230 + mrSges(1,1);
t152 = 2 * t175;
t148 = t159 * g(3);
t87 = t321 + 0.2e1 * (t162 + t378) * t220 - 0.2e1 * t138 + 0.2e1 * t151 + t241;
t86 = t322 + 0.2e1 * (t161 + t379) * t217 - 0.2e1 * t137 + 0.2e1 * t150 + t241;
t85 = t323 + 0.2e1 * (t160 + t380) * t214 - 0.2e1 * t136 + 0.2e1 * t149 + t241;
t63 = (t123 * t270 + t90 * t340) * t234;
t62 = (t122 * t271 + t89 * t342) * t234;
t61 = (t121 * t272 + t88 * t344) * t234;
t60 = (t90 * t270 + t87 * t340) * t234;
t59 = (t89 * t271 + t86 * t342) * t234;
t58 = (t88 * t272 + t85 * t344) * t234;
t57 = Ifges(3,3) * t261 + (t105 * t351 + t84 * t264) * t330;
t56 = Ifges(3,3) * t258 + (t104 * t351 + t83 * t264) * t330;
t55 = Ifges(3,3) * t262 + (t103 * t352 + t82 * t265) * t331;
t54 = Ifges(3,3) * t259 + (t102 * t352 + t81 * t265) * t331;
t53 = Ifges(3,3) * t263 + (t101 * t353 + t80 * t266) * t332;
t52 = Ifges(3,3) * t260 + (t100 * t353 + t79 * t266) * t332;
t48 = t144 * t261 + (t84 * t267 + t90 * t360) * t330;
t47 = t144 * t258 + (t83 * t267 + t90 * t361) * t330;
t46 = t143 * t262 + (t82 * t268 + t89 * t362) * t331;
t45 = t143 * t259 + (t81 * t268 + t89 * t363) * t331;
t44 = t142 * t263 + (t80 * t269 + t88 * t364) * t332;
t43 = t142 * t260 + (t79 * t269 + t88 * t365) * t332;
t42 = t117 * t261 + (t84 * t276 + t87 * t360) * t330;
t41 = t117 * t258 + (t83 * t276 + t87 * t361) * t330;
t40 = t116 * t262 + (t82 * t279 + t86 * t362) * t331;
t39 = t116 * t259 + (t81 * t279 + t86 * t363) * t331;
t38 = t115 * t263 + (t80 * t282 + t85 * t364) * t332;
t37 = t115 * t260 + (t79 * t282 + t85 * t365) * t332;
t30 = t78 / 0.2e1 + t77 / 0.2e1 + t99 / 0.2e1 + t69;
t29 = t76 / 0.2e1 + t75 / 0.2e1 + t98 / 0.2e1 + t68;
t28 = t74 / 0.2e1 + t73 / 0.2e1 + t97 / 0.2e1 + t67;
t15 = (-t408 + (-t111 * t192 - t220 * t33) * pkin(2)) * t330;
t14 = (-t409 + (-t110 * t189 - t217 * t32) * pkin(2)) * t331;
t13 = (-t410 + (-t109 * t186 - t214 * t31) * pkin(2)) * t332;
t12 = (-t231 * t367 + (pkin(1) * t69 + (0.2e1 * t30 * t308 - t273) * pkin(2)) * pkin(1)) * t192 * t69 * t252 + (-pkin(2) * t367 + (-t308 * t393 - t273) * pkin(1)) * t51 * t253 + ((pkin(1) * t311 * t393 + pkin(2) * t114) * t220 + pkin(1) * t354) * t193 * t114 * t330;
t11 = (-t231 * t369 + (pkin(1) * t68 + (0.2e1 * t29 * t309 - t274) * pkin(2)) * pkin(1)) * t189 * t68 * t254 + (-pkin(2) * t369 + (-t309 * t394 - t274) * pkin(1)) * t50 * t255 + ((pkin(1) * t312 * t394 + pkin(2) * t113) * t217 + pkin(1) * t355) * t190 * t113 * t331;
t10 = (-t231 * t371 + (pkin(1) * t67 + (0.2e1 * t28 * t310 - t275) * pkin(2)) * pkin(1)) * t186 * t67 * t256 + (-pkin(2) * t371 + (-t310 * t395 - t275) * pkin(1)) * t49 * t257 + ((pkin(1) * t313 * t395 + pkin(2) * t112) * t214 + pkin(1) * t356) * t187 * t112 * t332;
t9 = -t117 * t15 - t144 * t12 + Ifges(3,3) * t192 * t357 + (mrSges(3,2) * t408 + t33 * t318) * t220 + mrSges(3,1) * t66 * t294 - (g(1) * t178 + g(2) * t181) * t242 + (g(3) * t168 + t129 * t165) * (mrSges(3,1) * t211 + mrSges(3,2) * t220) + (-0.2e1 * t191 + 0.1e1) * t33 * Ifges(3,4);
t8 = -t116 * t14 - t143 * t11 + Ifges(3,3) * t189 * t358 + (mrSges(3,2) * t409 + t32 * t319) * t217 + mrSges(3,1) * t65 * t295 - (g(1) * t177 + g(2) * t180) * t243 + (g(3) * t167 + t128 * t164) * (mrSges(3,1) * t208 + mrSges(3,2) * t217) + (-0.2e1 * t188 + 0.1e1) * t32 * Ifges(3,4);
t7 = -t115 * t13 - t142 * t10 + Ifges(3,3) * t186 * t359 + (mrSges(3,2) * t410 + t31 * t320) * t214 + mrSges(3,1) * t64 * t296 - (g(1) * t176 + g(2) * t179) * t244 + (g(3) * t166 + t127 * t163) * (mrSges(3,1) * t205 + mrSges(3,2) * t214) + (-0.2e1 * t185 + 0.1e1) * t31 * Ifges(3,4);
t6 = -t90 * t15 - t123 * t12 - t175 * t168 + (t242 * g(3) + t226) * t165 + (-t198 * t165 - t235 * t168) * t129 + (t144 * t192 - Ifges(3,6)) * t357 + (t235 * t212 - t315) * t66 * pkin(1) + t238;
t5 = -t89 * t14 - t122 * t11 - t175 * t167 + (t243 * g(3) + t226) * t164 + (-t198 * t164 - t236 * t167) * t128 + (t143 * t189 - Ifges(3,6)) * t358 + (t236 * t209 - t316) * t65 * pkin(1) + t239;
t4 = -t88 * t13 - t121 * t10 - t175 * t166 + (t244 * g(3) + t226) * t163 + (-t198 * t163 - t237 * t166) * t127 + (t142 * t186 - Ifges(3,6)) * t359 + (t237 * t206 - t317) * t64 * pkin(1) + t240;
t3 = -t87 * t15 - t90 * t12 + (-t223 - t390) * t158 / 0.2e1 + (t225 - t384) * sin(t174) / 0.2e1 + (t223 - t390) * t157 / 0.2e1 + (t225 + t384) * sin(t173) / 0.2e1 + (t129 * t406 - t152) * t168 / 0.2e1 + t213 * (mrSges(1,2) * t129 + t148) + (-Ifges(3,6) + t351) * t357 + (-t129 * t198 + t226) * t165 + (-t129 * t159 + t224) * t222 + ((-mrSges(3,1) * t372 + mrSges(3,2) * t285) * t220 + (mrSges(3,1) * t285 + mrSges(3,2) * t372) * t211 + (-mrSges(2,1) * t212 + t315) * t51 * t30) * t411 + t238;
t2 = -t86 * t14 - t89 * t11 + (-t223 - t391) * t156 / 0.2e1 + (t225 - t385) * sin(t172) / 0.2e1 + (t223 - t391) * t155 / 0.2e1 + (t225 + t385) * sin(t171) / 0.2e1 + (t128 * t406 - t152) * t167 / 0.2e1 + t210 * (mrSges(1,2) * t128 + t148) + (-Ifges(3,6) + t352) * t358 + (-t128 * t198 + t226) * t164 + (-t128 * t159 + t224) * t219 + ((-mrSges(3,1) * t373 + mrSges(3,2) * t286) * t217 + (mrSges(3,1) * t286 + mrSges(3,2) * t373) * t208 + (-mrSges(2,1) * t209 + t316) * t50 * t29) * t411 + t239;
t1 = -t85 * t13 - t88 * t10 + (-t223 - t392) * t154 / 0.2e1 + (t225 - t386) * sin(t170) / 0.2e1 + (t223 - t392) * t153 / 0.2e1 + (t225 + t386) * sin(t169) / 0.2e1 + (t127 * t406 - t152) * t166 / 0.2e1 + t207 * (mrSges(1,2) * t127 + t148) + (-Ifges(3,6) + t353) * t359 + (-t127 * t198 + t226) * t163 + (-t127 * t159 + t224) * t216 + ((-mrSges(3,1) * t374 + mrSges(3,2) * t287) * t214 + (mrSges(3,1) * t287 + mrSges(3,2) * t374) * t205 + (-mrSges(2,1) * t206 + t317) * t49 * t28) * t411 + t240;
t16 = [(-g(1) + t204) * m(4) + ((t181 * t203 * t57 + (t204 * t57 + t9) * t178) * t192 + (t180 * t203 * t55 + (t204 * t55 + t8) * t177) * t189 + (t179 * t203 * t53 + (t204 * t53 + t7) * t176) * t186) * t232 + ((t48 * t348 + t46 * t349 + t44 * t350) * t314 + ((t48 * t277 + t42 * t360) * t204 + (t48 * t278 + t42 * t361) * t203 + t42 * t339 + t3 * t360 + t6 * t277) * t184 + ((t46 * t280 + t40 * t362) * t204 + (t46 * t281 + t40 * t363) * t203 + t40 * t341 + t2 * t362 + t5 * t280) * t183 + ((t44 * t283 + t38 * t364) * t204 + (t44 * t284 + t38 * t365) * t203 + t38 * t343 + t1 * t364 + t4 * t283) * t182) * t234; (-g(2) + t203) * m(4) + ((t178 * t204 * t56 + (t203 * t56 + t9) * t181) * t192 + (t177 * t204 * t54 + (t203 * t54 + t8) * t180) * t189 + (t176 * t204 * t52 + (t203 * t52 + t7) * t179) * t186) * t232 + ((t47 * t348 + t45 * t349 + t43 * t350) * t314 + ((t47 * t277 + t41 * t360) * t204 + (t47 * t278 + t41 * t361) * t203 + t41 * t339 + t3 * t361 + t6 * t278) * t184 + ((t45 * t280 + t39 * t362) * t204 + (t45 * t281 + t39 * t363) * t203 + t39 * t341 + t2 * t363 + t5 * t281) * t183 + ((t43 * t283 + t37 * t364) * t204 + (t43 * t284 + t37 * t365) * t203 + t37 * t343 + t1 * t365 + t4 * t284) * t182) * t234; (-g(3) + t202) * m(4) + ((t176 * t370 + t177 * t368 + t178 * t366) * t204 + (t179 * t370 + t180 * t368 + t181 * t366) * t203) * t232 + ((t4 * t350 + t5 * t349 + t6 * t348 + (t63 * t348 + t62 * t349 + t61 * t350) * t202) * t232 + ((t63 * t277 + t60 * t360) * t204 + (t63 * t278 + t60 * t361) * t203 + (t202 * t60 + t3) * t168) * t184 + ((t62 * t280 + t59 * t362) * t204 + (t62 * t281 + t59 * t363) * t203 + (t202 * t59 + t2) * t167) * t183 + ((t61 * t283 + t58 * t364) * t204 + (t61 * t284 + t58 * t365) * t203 + (t202 * t58 + t1) * t166) * t182) * t234;];
tauX  = t16;
