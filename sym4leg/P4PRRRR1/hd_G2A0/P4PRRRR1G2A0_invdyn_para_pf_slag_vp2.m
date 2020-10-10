% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
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
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:38
% EndTime: 2020-08-07 10:58:43
% DurationCPUTime: 5.11s
% Computational Cost: add. (6647->459), mult. (13057->849), div. (5940->13), fcn. (13506->26), ass. (0->364)
t192 = sin(qJ(3,4));
t223 = xP(4);
t172 = sin(t223);
t173 = cos(t223);
t226 = koppelP(4,2);
t230 = koppelP(4,1);
t135 = t172 * t230 + t173 * t226;
t139 = -t172 * t226 + t173 * t230;
t198 = legFrame(4,2);
t164 = sin(t198);
t168 = cos(t198);
t219 = xDP(4);
t221 = xDP(2);
t222 = xDP(1);
t282 = (t139 * t219 + t221) * t164 - (-t135 * t219 + t222) * t168;
t194 = cos(qJ(3,4));
t176 = 0.1e1 / t194;
t234 = 0.1e1 / pkin(2);
t328 = t176 * t234;
t81 = t282 * t328;
t77 = t81 ^ 2;
t407 = t192 * t77;
t206 = sin(qJ(3,3));
t227 = koppelP(3,2);
t231 = koppelP(3,1);
t136 = t172 * t231 + t173 * t227;
t140 = -t172 * t227 + t173 * t231;
t199 = legFrame(3,2);
t165 = sin(t199);
t169 = cos(t199);
t281 = (t140 * t219 + t221) * t165 - (-t136 * t219 + t222) * t169;
t212 = cos(qJ(3,3));
t183 = 0.1e1 / t212;
t321 = t183 * t234;
t82 = t281 * t321;
t78 = t82 ^ 2;
t406 = t206 * t78;
t208 = sin(qJ(3,2));
t228 = koppelP(2,2);
t232 = koppelP(2,1);
t137 = t172 * t232 + t173 * t228;
t141 = -t172 * t228 + t173 * t232;
t200 = legFrame(2,2);
t166 = sin(t200);
t170 = cos(t200);
t280 = (t141 * t219 + t221) * t166 - (-t137 * t219 + t222) * t170;
t214 = cos(qJ(3,2));
t186 = 0.1e1 / t214;
t320 = t186 * t234;
t83 = t280 * t320;
t79 = t83 ^ 2;
t405 = t208 * t79;
t210 = sin(qJ(3,1));
t229 = koppelP(1,2);
t233 = koppelP(1,1);
t138 = t172 * t233 + t173 * t229;
t142 = -t172 * t229 + t173 * t233;
t201 = legFrame(1,2);
t167 = sin(t201);
t171 = cos(t201);
t279 = (t142 * t219 + t221) * t167 - (-t138 * t219 + t222) * t171;
t216 = cos(qJ(3,1));
t189 = 0.1e1 / t216;
t319 = t189 * t234;
t84 = t279 * t319;
t80 = t84 ^ 2;
t404 = t210 * t80;
t191 = t219 ^ 2;
t202 = xDDP(4);
t205 = xDDP(1);
t101 = -t138 * t202 - t142 * t191 + t205;
t351 = t101 * t171;
t204 = xDDP(2);
t97 = -t138 * t191 + t142 * t202 + t204;
t377 = t167 * t97;
t403 = t351 - t377;
t100 = -t137 * t202 - t141 * t191 + t205;
t353 = t100 * t170;
t96 = -t137 * t191 + t141 * t202 + t204;
t378 = t166 * t96;
t402 = t353 - t378;
t99 = -t136 * t202 - t140 * t191 + t205;
t375 = t169 * t99;
t95 = -t136 * t191 + t140 * t202 + t204;
t379 = t165 * t95;
t401 = t375 - t379;
t98 = -t135 * t202 - t139 * t191 + t205;
t376 = t168 * t98;
t94 = -t135 * t191 + t139 * t202 + t204;
t380 = t164 * t94;
t400 = t376 - t380;
t188 = t216 ^ 2;
t399 = (0.2e1 * t188 - 0.1e1) * Ifges(3,4);
t185 = t214 ^ 2;
t398 = (0.2e1 * t185 - 0.1e1) * Ifges(3,4);
t182 = t212 ^ 2;
t397 = (0.2e1 * t182 - 0.1e1) * Ifges(3,4);
t175 = t194 ^ 2;
t396 = (0.2e1 * t175 - 0.1e1) * Ifges(3,4);
t395 = 0.2e1 * Ifges(3,4);
t390 = Ifges(3,5) / 0.2e1;
t389 = -Ifges(3,6) / 0.2e1;
t197 = mrSges(2,2) - mrSges(3,3);
t388 = t197 / 0.2e1;
t387 = Ifges(3,1) + Ifges(2,3);
t193 = sin(qJ(2,4));
t318 = t193 * t194;
t121 = t164 * t318 - t168 * t192;
t386 = t121 * t98;
t122 = t164 * t192 + t168 * t318;
t385 = t122 * t94;
t207 = sin(qJ(2,3));
t313 = t207 * t212;
t123 = t165 * t313 - t169 * t206;
t384 = t123 * t99;
t126 = t165 * t206 + t169 * t313;
t383 = t126 * t95;
t209 = sin(qJ(2,2));
t312 = t209 * t214;
t127 = t166 * t208 + t170 * t312;
t382 = t127 * t96;
t211 = sin(qJ(2,1));
t311 = t211 * t216;
t128 = t167 * t210 + t171 * t311;
t381 = t128 * t97;
t374 = t176 * t77;
t373 = t183 * t78;
t372 = t186 * t79;
t371 = t189 * t80;
t174 = 0.1e1 / t193;
t220 = xDP(3);
t177 = 0.1e1 / t194 ^ 2;
t195 = cos(qJ(2,4));
t301 = t177 * t192 * t195;
t54 = (-t176 * t220 - t282 * t301) * t234 * t174;
t370 = t192 * t54;
t369 = t192 * t81;
t53 = t54 ^ 2;
t368 = t194 * t53;
t367 = t195 * t54;
t179 = 0.1e1 / t207;
t184 = 0.1e1 / t212 ^ 2;
t213 = cos(qJ(2,3));
t297 = t184 * t206 * t213;
t58 = (-t183 * t220 - t281 * t297) * t234 * t179;
t366 = t206 * t58;
t365 = t206 * t82;
t180 = 0.1e1 / t209;
t187 = 0.1e1 / t214 ^ 2;
t215 = cos(qJ(2,2));
t296 = t187 * t208 * t215;
t59 = (-t186 * t220 - t280 * t296) * t234 * t180;
t364 = t208 * t59;
t363 = t208 * t83;
t181 = 0.1e1 / t211;
t190 = 0.1e1 / t216 ^ 2;
t217 = cos(qJ(2,1));
t295 = t190 * t210 * t217;
t60 = (-t189 * t220 - t279 * t295) * t234 * t181;
t362 = t210 * t60;
t361 = t210 * t84;
t55 = t58 ^ 2;
t360 = t212 * t55;
t359 = t213 * t58;
t56 = t59 ^ 2;
t358 = t214 * t56;
t357 = t215 * t59;
t57 = t60 ^ 2;
t356 = t216 * t57;
t355 = t217 * t60;
t124 = t166 * t312 - t170 * t208;
t354 = t100 * t124;
t125 = t167 * t311 - t171 * t210;
t352 = t101 * t125;
t151 = Ifges(3,5) * t192 + Ifges(3,6) * t194;
t350 = t151 * t176;
t152 = mrSges(3,1) * t192 + mrSges(3,2) * t194;
t349 = t152 * t176;
t348 = t152 * t193;
t153 = Ifges(3,5) * t206 + Ifges(3,6) * t212;
t347 = t153 * t183;
t154 = Ifges(3,5) * t208 + Ifges(3,6) * t214;
t346 = t154 * t186;
t155 = Ifges(3,5) * t210 + Ifges(3,6) * t216;
t345 = t155 * t189;
t156 = mrSges(3,1) * t206 + mrSges(3,2) * t212;
t344 = t156 * t183;
t343 = t156 * t207;
t157 = mrSges(3,1) * t208 + mrSges(3,2) * t214;
t342 = t157 * t186;
t341 = t157 * t209;
t158 = mrSges(3,1) * t210 + mrSges(3,2) * t216;
t340 = t158 * t189;
t339 = t158 * t211;
t338 = t164 * t234;
t337 = t165 * t234;
t336 = t166 * t234;
t335 = t167 * t234;
t334 = t168 * t234;
t333 = t169 * t234;
t332 = t170 * t234;
t331 = t171 * t234;
t330 = t174 * t176;
t329 = t174 * t195;
t327 = t179 * t183;
t326 = t179 * t213;
t325 = t180 * t186;
t324 = t180 * t215;
t323 = t181 * t189;
t322 = t181 * t217;
t203 = xDDP(3);
t317 = t195 * t203;
t316 = t203 * t213;
t315 = t203 * t215;
t314 = t203 * t217;
t290 = mrSges(3,1) * t194 - mrSges(3,2) * t192;
t274 = mrSges(2,1) + t290;
t109 = -t193 * t197 + t274 * t195;
t310 = t109 * t330;
t289 = mrSges(3,1) * t212 - mrSges(3,2) * t206;
t273 = mrSges(2,1) + t289;
t110 = -t207 * t197 + t273 * t213;
t309 = t110 * t327;
t288 = mrSges(3,1) * t214 - mrSges(3,2) * t208;
t272 = mrSges(2,1) + t288;
t111 = -t209 * t197 + t272 * t215;
t308 = t111 * t325;
t287 = mrSges(3,1) * t216 - mrSges(3,2) * t210;
t271 = mrSges(2,1) + t287;
t112 = -t211 * t197 + t271 * t217;
t307 = t112 * t323;
t306 = t176 * t348;
t305 = t183 * t343;
t304 = t186 * t341;
t303 = t189 * t339;
t178 = m(1) + m(2) + m(3);
t302 = t178 * t330;
t300 = t178 * t327;
t299 = t178 * t325;
t298 = t178 * t323;
t294 = t174 * t301;
t293 = t179 * t297;
t292 = t180 * t296;
t291 = t181 * t295;
t143 = g(1) * t164 + g(2) * t168;
t286 = g(3) * t195 + t143 * t193;
t144 = g(1) * t165 + g(2) * t169;
t285 = g(3) * t213 + t144 * t207;
t145 = g(1) * t166 + g(2) * t170;
t284 = g(3) * t215 + t145 * t209;
t146 = g(1) * t167 + g(2) * t171;
t283 = g(3) * t217 + t146 * t211;
t278 = t135 * t168 + t139 * t164;
t277 = t136 * t169 + t140 * t165;
t276 = t137 * t170 + t141 * t166;
t275 = t138 * t171 + t142 * t167;
t270 = t176 * (t385 + t386);
t269 = t183 * (t383 + t384);
t268 = t186 * (t354 + t382);
t267 = t189 * (t352 + t381);
t196 = Ifges(3,1) - Ifges(3,2);
t129 = t192 * t194 * t395 - t175 * t196 + t387;
t13 = ((-t193 * t369 + t194 * t367) * t176 * t54 + (t195 * t81 - t318 * t370) * t177 * t81) * t174;
t163 = t197 * g(3);
t218 = mrSges(2,1) * g(3);
t29 = (-t368 - t374) * t174 * pkin(2);
t5 = -t109 * t29 - t129 * t13 + t350 * t407 + 0.2e1 * t81 * ((t196 * t370 + t81 * t390) * t194 + t369 * t389 + t54 * t396) + (-t274 * t143 + t163) * t195 + t193 * (t290 * g(3) + t143 * t197 + t218);
t147 = g(1) * t168 - g(2) * t164;
t9 = t29 * t348 - t151 * t13 + (mrSges(3,1) * t147 + t286 * mrSges(3,2)) * t194 + (t286 * mrSges(3,1) - mrSges(3,2) * t147 + Ifges(3,3) * t374 - t196 * t368) * t192 - t53 * t396;
t266 = t176 * t9 - t5 * t294;
t14 = ((-t207 * t365 + t212 * t359) * t183 * t58 + (t213 * t82 - t313 * t366) * t184 * t82) * t179;
t148 = g(1) * t169 - g(2) * t165;
t30 = (-t360 - t373) * t179 * pkin(2);
t10 = t30 * t343 - t153 * t14 + (mrSges(3,1) * t148 + t285 * mrSges(3,2)) * t212 + (t285 * mrSges(3,1) - mrSges(3,2) * t148 + Ifges(3,3) * t373 - t196 * t360) * t206 - t55 * t397;
t130 = t206 * t212 * t395 - t182 * t196 + t387;
t6 = -t110 * t30 - t130 * t14 + t347 * t406 + 0.2e1 * t82 * ((t196 * t366 + t82 * t390) * t212 + t365 * t389 + t58 * t397) + (-t273 * t144 + t163) * t213 + t207 * (t289 * g(3) + t144 * t197 + t218);
t265 = t183 * t10 - t6 * t293;
t149 = g(1) * t170 - g(2) * t166;
t15 = ((-t209 * t363 + t214 * t357) * t186 * t59 + (t215 * t83 - t312 * t364) * t187 * t83) * t180;
t31 = (-t358 - t372) * t180 * pkin(2);
t11 = t31 * t341 - t154 * t15 + (mrSges(3,1) * t149 + t284 * mrSges(3,2)) * t214 + (t284 * mrSges(3,1) - mrSges(3,2) * t149 + Ifges(3,3) * t372 - t196 * t358) * t208 - t56 * t398;
t131 = t208 * t214 * t395 - t185 * t196 + t387;
t7 = -t111 * t31 - t131 * t15 + t346 * t405 + 0.2e1 * t83 * ((t196 * t364 + t83 * t390) * t214 + t363 * t389 + t59 * t398) + (-t272 * t145 + t163) * t215 + t209 * (t288 * g(3) + t145 * t197 + t218);
t264 = t186 * t11 - t7 * t292;
t150 = g(1) * t171 - g(2) * t167;
t16 = ((-t211 * t361 + t216 * t355) * t189 * t60 + (t217 * t84 - t311 * t362) * t190 * t84) * t181;
t32 = (-t356 - t371) * t181 * pkin(2);
t12 = t32 * t339 - t155 * t16 + (mrSges(3,1) * t150 + t283 * mrSges(3,2)) * t216 + (t283 * mrSges(3,1) - mrSges(3,2) * t150 + Ifges(3,3) * t371 - t196 * t356) * t210 - t57 * t399;
t132 = t210 * t216 * t395 - t188 * t196 + t387;
t8 = -t112 * t32 - t132 * t16 + t345 * t404 + 0.2e1 * t84 * ((t196 * t362 + t84 * t390) * t216 + t361 * t389 + t60 * t399) + (-t271 * t146 + t163) * t217 + t211 * (t287 * g(3) + t146 * t197 + t218);
t263 = t189 * t12 - t8 * t291;
t250 = -Ifges(3,3) * t176 + t151 * t294;
t242 = t129 * t294 - t350;
t37 = t121 * t310 + t242 * t334;
t261 = -t176 * (-t121 * t349 + t250 * t334) + t37 * t294;
t38 = t122 * t310 - t242 * t338;
t260 = -t176 * (-t122 * t349 - t250 * t338) + t38 * t294;
t249 = -Ifges(3,3) * t183 + t153 * t293;
t241 = t130 * t293 - t347;
t39 = t123 * t309 + t241 * t333;
t258 = -t183 * (-t123 * t344 + t249 * t333) + t39 * t293;
t42 = t126 * t309 - t241 * t337;
t257 = -t183 * (-t126 * t344 - t249 * t337) + t42 * t293;
t248 = -Ifges(3,3) * t186 + t154 * t292;
t240 = t131 * t292 - t346;
t40 = t124 * t308 + t240 * t332;
t255 = -t186 * (-t124 * t342 + t248 * t332) + t40 * t292;
t43 = t127 * t308 - t240 * t336;
t254 = -t186 * (-t127 * t342 - t248 * t336) + t43 * t292;
t247 = -Ifges(3,3) * t189 + t155 * t291;
t239 = t132 * t291 - t345;
t41 = t125 * t307 + t239 * t331;
t252 = -t189 * (-t125 * t340 + t247 * t331) + t41 * t291;
t44 = t128 * t307 - t239 * t335;
t251 = -t189 * (-t128 * t340 - t247 * t335) + t44 * t291;
t238 = t109 * t294 + t306;
t237 = t110 * t293 + t305;
t236 = t111 * t292 + t304;
t235 = t112 * t291 + t303;
t225 = mrSges(4,1);
t224 = mrSges(4,2);
t134 = -t172 * t224 + t173 * t225;
t133 = t172 * t225 + t173 * t224;
t104 = (-t112 * t319 + t178 * t217) * t181;
t103 = (-t111 * t320 + t178 * t215) * t180;
t102 = (-t110 * t321 + t178 * t213) * t179;
t93 = (-t109 * t328 + t178 * t195) * t174;
t92 = t275 * t319;
t91 = t276 * t320;
t90 = t277 * t321;
t89 = t278 * t328;
t88 = (t112 * t217 - t132 * t319) * t181;
t87 = (t111 * t215 - t131 * t320) * t180;
t86 = (t110 * t213 - t130 * t321) * t179;
t85 = (t109 * t195 - t129 * t328) * t174;
t76 = t275 * t234 * t291;
t75 = t276 * t234 * t292;
t74 = t277 * t234 * t293;
t73 = t278 * t234 * t294;
t64 = (-t125 * t138 + t128 * t142) * t323;
t63 = (-t124 * t137 + t127 * t141) * t325;
t62 = (-t123 * t136 + t126 * t140) * t327;
t61 = (-t121 * t135 + t122 * t139) * t330;
t52 = t128 * t298 - t235 * t335;
t51 = t127 * t299 - t236 * t336;
t50 = t126 * t300 - t237 * t337;
t49 = t125 * t298 + t235 * t331;
t48 = t124 * t299 + t236 * t332;
t47 = t123 * t300 + t237 * t333;
t46 = t122 * t302 - t238 * t338;
t45 = t121 * t302 + t238 * t334;
t24 = -t112 * t76 + t178 * t64 - t92 * t339;
t23 = -t111 * t75 + t178 * t63 - t91 * t341;
t22 = -t110 * t74 + t178 * t62 - t90 * t343;
t21 = -t109 * t73 + t178 * t61 - t89 * t348;
t20 = t112 * t64 - t132 * t76 + t155 * t92;
t19 = t111 * t63 - t131 * t75 + t154 * t91;
t18 = t110 * t62 - t130 * t74 + t153 * t90;
t17 = t109 * t61 - t129 * t73 + t151 * t89;
t4 = -t112 * t16 - t303 * t404 + (-mrSges(2,1) * t57 - t287 * (t57 + t80)) * t211 - 0.2e1 * (t158 * t84 + t60 * t388) * t355 + (-t32 - t146) * t178;
t3 = -t111 * t15 - t304 * t405 + (-mrSges(2,1) * t56 - t288 * (t56 + t79)) * t209 - 0.2e1 * (t157 * t83 + t59 * t388) * t357 + (-t31 - t145) * t178;
t2 = -t110 * t14 - t305 * t406 + (-mrSges(2,1) * t55 - t289 * (t55 + t78)) * t207 - 0.2e1 * (t156 * t82 + t58 * t388) * t359 + (-t30 - t144) * t178;
t1 = -t109 * t13 - t306 * t407 + (-mrSges(2,1) * t53 - t290 * (t53 + t77)) * t193 - 0.2e1 * (t152 * t81 + t54 * t388) * t367 + (-t29 - t143) * t178;
t25 = [-t133 * t202 - t191 * t134 + (t205 - g(1)) * m(4) + (t49 * t314 + (t49 * t381 + (t101 * t49 + t4) * t125) * t189) * t181 + (t48 * t315 + (t48 * t382 + (t100 * t48 + t3) * t124) * t186) * t180 + (t47 * t316 + (t47 * t383 + (t47 * t99 + t2) * t123) * t183) * t179 + (t45 * t317 + (t45 * t385 + (t45 * t98 + t1) * t121) * t176) * t174 + (-t252 * t377 - t255 * t378 - t258 * t379 - t261 * t380 + (-t41 * t323 - t40 * t325 - t39 * t327 - t37 * t330) * t203 + (t101 * t252 - t263) * t171 + (t100 * t255 - t264) * t170 + (t258 * t99 - t265) * t169 + (t261 * t98 - t266) * t168) * t234; -t191 * t133 + t134 * t202 + (t204 - g(2)) * m(4) + (t52 * t314 + (t52 * t352 + (t52 * t97 + t4) * t128) * t189) * t181 + (t51 * t315 + (t51 * t354 + (t51 * t96 + t3) * t127) * t186) * t180 + (t50 * t316 + (t50 * t384 + (t50 * t95 + t2) * t126) * t183) * t179 + (t46 * t317 + (t46 * t386 + (t46 * t94 + t1) * t122) * t176) * t174 + (t251 * t351 + t254 * t353 + t257 * t375 + t260 * t376 + (-t44 * t323 - t43 * t325 - t42 * t327 - t38 * t330) * t203 + (-t251 * t97 + t263) * t167 + (-t254 * t96 + t264) * t166 + (-t257 * t95 + t265) * t165 + (-t260 * t94 + t266) * t164) * t234; -g(3) * m(4) + (t104 * t267 + t217 * t4) * t181 + (t103 * t268 + t215 * t3) * t180 + (t102 * t269 + t213 * t2) * t179 + (t195 * t1 + t270 * t93) * t174 + (t102 * t326 + t103 * t324 + t104 * t322 + t93 * t329 + m(4)) * t203 + (-t5 * t330 - t6 * t327 - t7 * t325 - t8 * t323 + (-t88 * t323 - t87 * t325 - t86 * t327 - t85 * t330) * t203 + t400 * (-(-t151 * t174 * t328 - t195 * t152) * t176 + t85 * t294) + t401 * (-(-t153 * t179 * t321 - t213 * t156) * t183 + t86 * t293) + t402 * (-(-t154 * t180 * t320 - t215 * t157) * t186 + t87 * t292) + t403 * (-(-t155 * t181 * t319 - t217 * t158) * t189 + t88 * t291)) * t234; -(-g(1) * t225 - g(2) * t224) * t172 + t173 * (g(1) * t224 - g(2) * t225) + Ifges(4,3) * t202 + t134 * t204 - t133 * t205 + t61 * t1 + t62 * t2 + t63 * t3 + t64 * t4 - t73 * t5 - t74 * t6 - t75 * t7 - t76 * t8 + t89 * t9 + t90 * t10 + t91 * t11 + t92 * t12 + t24 * t181 * t267 + t23 * t180 * t268 + t22 * t179 * t269 + t21 * t174 * t270 + (t21 * t329 + t22 * t326 + t23 * t324 + t24 * t322) * t203 + ((-t17 * t330 - t18 * t327 - t19 * t325 - t20 * t323) * t203 + t400 * (t17 * t294 - t176 * (Ifges(3,3) * t89 - t151 * t73 - t61 * t348)) + t401 * (t18 * t293 - t183 * (Ifges(3,3) * t90 - t153 * t74 - t62 * t343)) + t402 * (-t186 * (Ifges(3,3) * t91 - t154 * t75 - t63 * t341) + t19 * t292) + t403 * (-t189 * (Ifges(3,3) * t92 - t155 * t76 - t64 * t339) + t20 * t291)) * t234;];
tauX  = t25;
