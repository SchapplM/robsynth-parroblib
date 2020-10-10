% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V2G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:26:47
% EndTime: 2020-08-06 19:26:58
% DurationCPUTime: 11.16s
% Computational Cost: add. (86913->640), mult. (102525->1045), div. (11898->6), fcn. (76671->18), ass. (0->381)
t209 = sin(qJ(2,3));
t434 = 0.2e1 * t209;
t211 = sin(qJ(2,2));
t433 = 0.2e1 * t211;
t213 = sin(qJ(2,1));
t432 = 0.2e1 * t213;
t399 = pkin(1) * t209;
t161 = qJ(3,3) + t399;
t169 = t209 * qJ(3,3);
t210 = sin(qJ(1,3));
t215 = cos(qJ(2,3));
t194 = t215 ^ 2;
t229 = pkin(2) + pkin(3);
t338 = (qJ(3,3) + t229) * (-qJ(3,3) + t229);
t276 = t194 * t338;
t313 = t229 * t215;
t216 = cos(qJ(1,3));
t228 = pkin(5) - pkin(6);
t317 = t228 * t216;
t103 = -t210 * t276 - ((0.2e1 * t169 + pkin(1)) * t210 - t317) * t313 - qJ(3,3) * (t161 * t210 - t209 * t317);
t225 = xDP(3);
t226 = xDP(2);
t227 = xDP(1);
t158 = t169 + pkin(1);
t430 = t158 + t313;
t137 = 0.1e1 / t430;
t231 = 0.1e1 / qJ(3,3);
t362 = t137 * t231;
t308 = pkin(1) * t216 + t210 * t228;
t121 = qJ(3,3) * t216 + t308 * t209;
t280 = t216 * t169;
t124 = 0.2e1 * t280 + t308;
t275 = t209 * t338;
t403 = pkin(1) * qJ(3,3);
t131 = -t275 + t403;
t203 = legFrame(3,2);
t172 = sin(t203);
t175 = cos(t203);
t274 = t216 * t338;
t416 = -0.2e1 * t229;
t295 = qJ(3,3) * t416;
t345 = t175 * t229;
t354 = t172 * t229;
t376 = qJ(3,3) * t172;
t76 = (-t172 * t274 + t175 * t295) * t194 + (-t124 * t354 - t131 * t175) * t215 - t121 * t376 + t161 * t345;
t375 = qJ(3,3) * t175;
t77 = (t172 * t295 + t175 * t274) * t194 + (t124 * t345 - t131 * t172) * t215 + t121 * t375 + t161 * t354;
t34 = (t103 * t225 + t226 * t76 + t227 * t77) * t362;
t109 = t430 * t210 - t317;
t365 = t109 * t215;
t125 = t280 + t308;
t324 = t216 * t229;
t333 = t209 * t229;
t88 = (-t172 * t324 - t375) * t194 + (-t125 * t172 + t175 * t333) * t215 + t175 * t161;
t89 = (t175 * t324 - t376) * t194 + (t125 * t175 + t172 * t333) * t215 + t172 * t161;
t61 = (-t225 * t365 + t226 * t88 + t227 * t89) * t362;
t386 = t229 * t61;
t25 = t34 - t386;
t398 = pkin(1) * t211;
t162 = qJ(3,2) + t398;
t170 = t211 * qJ(3,2);
t212 = sin(qJ(1,2));
t217 = cos(qJ(2,2));
t195 = t217 ^ 2;
t337 = (qJ(3,2) + t229) * (-qJ(3,2) + t229);
t273 = t195 * t337;
t312 = t229 * t217;
t218 = cos(qJ(1,2));
t316 = t228 * t218;
t104 = -t212 * t273 - ((0.2e1 * t170 + pkin(1)) * t212 - t316) * t312 - qJ(3,2) * (t162 * t212 - t211 * t316);
t159 = t170 + pkin(1);
t429 = t159 + t312;
t138 = 0.1e1 / t429;
t233 = 0.1e1 / qJ(3,2);
t361 = t138 * t233;
t307 = pkin(1) * t218 + t212 * t228;
t122 = qJ(3,2) * t218 + t307 * t211;
t281 = t218 * t170;
t126 = 0.2e1 * t281 + t307;
t272 = t211 * t337;
t404 = pkin(1) * qJ(3,2);
t132 = -t272 + t404;
t204 = legFrame(2,2);
t173 = sin(t204);
t176 = cos(t204);
t271 = t218 * t337;
t296 = qJ(3,2) * t416;
t342 = t176 * t229;
t351 = t173 * t229;
t378 = qJ(3,2) * t173;
t78 = (-t173 * t271 + t176 * t296) * t195 + (-t126 * t351 - t132 * t176) * t217 - t122 * t378 + t162 * t342;
t377 = qJ(3,2) * t176;
t79 = (t173 * t296 + t176 * t271) * t195 + (t126 * t342 - t132 * t173) * t217 + t122 * t377 + t162 * t351;
t35 = (t104 * t225 + t226 * t78 + t227 * t79) * t361;
t110 = t429 * t212 - t316;
t364 = t110 * t217;
t127 = t281 + t307;
t321 = t218 * t229;
t330 = t211 * t229;
t90 = (-t173 * t321 - t377) * t195 + (-t127 * t173 + t176 * t330) * t217 + t176 * t162;
t91 = (t176 * t321 - t378) * t195 + (t127 * t176 + t173 * t330) * t217 + t173 * t162;
t62 = (-t225 * t364 + t226 * t90 + t227 * t91) * t361;
t385 = t229 * t62;
t26 = t35 - t385;
t397 = pkin(1) * t213;
t163 = qJ(3,1) + t397;
t171 = t213 * qJ(3,1);
t214 = sin(qJ(1,1));
t219 = cos(qJ(2,1));
t196 = t219 ^ 2;
t336 = (qJ(3,1) + t229) * (-qJ(3,1) + t229);
t270 = t196 * t336;
t311 = t229 * t219;
t220 = cos(qJ(1,1));
t315 = t228 * t220;
t105 = -t214 * t270 - ((0.2e1 * t171 + pkin(1)) * t214 - t315) * t311 - qJ(3,1) * (t163 * t214 - t213 * t315);
t160 = t171 + pkin(1);
t428 = t160 + t311;
t139 = 0.1e1 / t428;
t235 = 0.1e1 / qJ(3,1);
t360 = t139 * t235;
t306 = pkin(1) * t220 + t214 * t228;
t123 = qJ(3,1) * t220 + t306 * t213;
t282 = t220 * t171;
t128 = 0.2e1 * t282 + t306;
t269 = t213 * t336;
t405 = pkin(1) * qJ(3,1);
t133 = -t269 + t405;
t205 = legFrame(1,2);
t174 = sin(t205);
t177 = cos(t205);
t268 = t220 * t336;
t297 = qJ(3,1) * t416;
t339 = t177 * t229;
t348 = t174 * t229;
t380 = qJ(3,1) * t174;
t80 = (-t174 * t268 + t177 * t297) * t196 + (-t128 * t348 - t133 * t177) * t219 - t123 * t380 + t163 * t339;
t379 = qJ(3,1) * t177;
t81 = (t174 * t297 + t177 * t268) * t196 + (t128 * t339 - t133 * t174) * t219 + t123 * t379 + t163 * t348;
t36 = (t105 * t225 + t226 * t80 + t227 * t81) * t360;
t111 = t428 * t214 - t315;
t363 = t111 * t219;
t129 = t282 + t306;
t318 = t220 * t229;
t327 = t213 * t229;
t92 = (-t174 * t318 - t379) * t196 + (-t129 * t174 + t177 * t327) * t219 + t177 * t163;
t93 = (t177 * t318 - t380) * t196 + (t129 * t177 + t174 * t327) * t219 + t174 * t163;
t63 = (-t225 * t363 + t226 * t92 + t227 * t93) * t360;
t384 = t229 * t63;
t27 = t36 - t384;
t305 = 0.2e1 * pkin(1);
t431 = 0.2e1 * t229;
t230 = qJ(3,3) ^ 2;
t423 = 2 * rSges(3,3);
t427 = qJ(3,3) * t423 + t230;
t232 = qJ(3,2) ^ 2;
t426 = qJ(3,2) * t423 + t232;
t234 = qJ(3,1) ^ 2;
t425 = qJ(3,1) * t423 + t234;
t424 = -0.2e1 * pkin(1);
t100 = (-t216 * t225 + (t172 * t226 - t175 * t227) * t210) * t137;
t422 = 0.2e1 * t100;
t101 = (-t218 * t225 + (t173 * t226 - t176 * t227) * t212) * t138;
t421 = 0.2e1 * t101;
t102 = (-t220 * t225 + (t174 * t226 - t177 * t227) * t214) * t139;
t420 = 0.2e1 * t102;
t419 = -0.2e1 * t194;
t418 = -0.2e1 * t195;
t417 = -0.2e1 * t196;
t415 = m(2) * rSges(2,1);
t224 = m(2) * rSges(2,2);
t223 = rSges(2,3) + pkin(5);
t222 = pkin(2) + rSges(3,1);
t414 = t194 - 0.1e1;
t413 = t195 - 0.1e1;
t412 = t196 - 0.1e1;
t411 = m(2) * t223;
t410 = m(3) * t137;
t409 = m(3) * t138;
t408 = m(3) * t139;
t221 = pkin(5) + rSges(3,2);
t407 = m(3) * t221;
t406 = m(3) * t222;
t200 = -qJ(3,3) - rSges(3,3);
t153 = m(3) * t200 + t224;
t402 = pkin(1) * t153;
t201 = -qJ(3,2) - rSges(3,3);
t154 = m(3) * t201 + t224;
t401 = pkin(1) * t154;
t202 = -qJ(3,1) - rSges(3,3);
t155 = m(3) * t202 + t224;
t400 = pkin(1) * t155;
t146 = m(1) * rSges(1,2) - t407 - t411;
t396 = g(3) * t146;
t22 = t228 * t25;
t286 = t100 * t403;
t395 = 0.2e1 * t286 + t22;
t23 = t228 * t26;
t287 = t101 * t404;
t394 = 0.2e1 * t287 + t23;
t24 = t228 * t27;
t288 = t102 * t405;
t393 = 0.2e1 * t288 + t24;
t392 = t200 * t34;
t391 = t201 * t35;
t390 = t202 * t36;
t334 = t209 * t228;
t279 = t100 * t334;
t389 = t215 * (t279 - t386);
t331 = t211 * t228;
t278 = t101 * t331;
t388 = t217 * (t278 - t385);
t328 = t213 * t228;
t277 = t102 * t328;
t387 = t219 * (t277 - t384);
t383 = t230 * t61;
t382 = t232 * t62;
t381 = t234 * t63;
t374 = t100 * t194;
t373 = t100 * t209;
t372 = t100 * t228;
t371 = t101 * t195;
t370 = t101 * t211;
t369 = t101 * t228;
t368 = t102 * t196;
t367 = t102 * t213;
t366 = t102 * t228;
t359 = t153 * t209;
t358 = t154 * t211;
t357 = t155 * t213;
t207 = xDDP(2);
t356 = t172 * t207;
t355 = t172 * t210;
t353 = t173 * t207;
t352 = t173 * t212;
t350 = t174 * t207;
t349 = t174 * t214;
t208 = xDDP(1);
t347 = t175 * t208;
t346 = t175 * t210;
t344 = t176 * t208;
t343 = t176 * t212;
t341 = t177 * t208;
t340 = t177 * t214;
t335 = t209 * t221;
t332 = t211 * t221;
t329 = t213 * t221;
t326 = t215 * t222;
t206 = xDDP(3);
t325 = t216 * t206;
t323 = t217 * t222;
t322 = t218 * t206;
t320 = t219 * t222;
t319 = t220 * t206;
t314 = t228 * t229;
t310 = -Icges(2,1) - Icges(3,1);
t309 = Icges(2,6) - Icges(3,6);
t304 = -0.2e1 * t407;
t303 = rSges(3,3) - t222;
t302 = rSges(3,3) + t222;
t300 = m(3) * t392;
t299 = m(3) * t391;
t298 = m(3) * t390;
t294 = t200 * t407;
t293 = t201 * t407;
t292 = t202 * t407;
t291 = m(3) * t335;
t290 = m(3) * t332;
t289 = m(3) * t329;
t285 = t215 * qJ(3,3) * t61;
t284 = t217 * qJ(3,2) * t62;
t283 = t63 * t219 * qJ(3,1);
t267 = t210 * t335;
t266 = t212 * t332;
t265 = t214 * t329;
t264 = -rSges(2,1) * t224 + Icges(2,4) - Icges(3,5);
t263 = rSges(2,2) * t411 - t309;
t236 = rSges(3,3) ^ 2;
t243 = pkin(1) ^ 2;
t262 = t221 ^ 2 + t236 + t243;
t260 = pkin(2) ^ 2 + t236 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t259 = (pkin(6) ^ 2) + t243 + ((-2 * pkin(6) + pkin(5)) * pkin(5));
t237 = rSges(2,2) ^ 2;
t239 = rSges(2,1) ^ 2;
t258 = (t237 + t239) * m(2) + Icges(3,2) + Icges(2,3);
t143 = g(1) * t175 - g(2) * t172;
t257 = -g(3) * t210 + t143 * t216;
t144 = g(1) * t176 - g(2) * t173;
t256 = -g(3) * t212 + t144 * t218;
t145 = g(1) * t177 - g(2) * t174;
t255 = -g(3) * t214 + t145 * t220;
t19 = t229 * t25 - t383;
t20 = t229 * t26 - t382;
t21 = t229 * t27 - t381;
t254 = Icges(2,2) + Icges(3,3) + (-t237 + t239) * m(2) + t310;
t156 = t406 + t415;
t253 = t156 * t215 - t359;
t252 = t156 * t217 - t358;
t251 = t156 * t219 - t357;
t250 = t221 * t406 - Icges(3,4) - Icges(2,5);
t249 = -t223 * t224 + t309;
t248 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (t223 ^ 2 + t237 + t243) * m(2) + Icges(1,3) - t310;
t247 = t223 * t415 + t250;
t199 = t229 ^ 2;
t152 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t151 = g(3) * t152;
t150 = t156 * t305;
t149 = t234 + t259;
t148 = t232 + t259;
t147 = t230 + t259;
t142 = g(1) * t174 + g(2) * t177;
t141 = g(1) * t173 + g(2) * t176;
t140 = g(1) * t172 + g(2) * t175;
t136 = -t202 * t406 + t264;
t135 = -t201 * t406 + t264;
t134 = -t200 * t406 + t264;
t130 = rSges(2,1) * t411 + t250;
t117 = -(qJ(3,1) + t302) * (qJ(3,1) + t303) * m(3) + t254;
t116 = -(qJ(3,2) + t302) * (qJ(3,2) + t303) * m(3) + t254;
t115 = -(qJ(3,3) + t302) * (qJ(3,3) + t303) * m(3) + t254;
t114 = (t260 + t425) * m(3) + t258;
t113 = (t260 + t426) * m(3) + t258;
t112 = (t260 + t427) * m(3) + t258;
t108 = (t249 - t292) * t219 - t213 * t247;
t107 = (t249 - t293) * t217 - t211 * t247;
t106 = (t249 - t294) * t215 - t209 * t247;
t99 = t102 ^ 2;
t98 = t101 ^ 2;
t97 = t100 ^ 2;
t84 = t117 * t196 + (t136 * t432 + t150) * t219 + t357 * t424 + (t262 + t425) * m(3) + t248;
t83 = t116 * t195 + (t135 * t433 + t150) * t217 + t358 * t424 + (t262 + t426) * m(3) + t248;
t82 = t115 * t194 + (t134 * t434 + t150) * t215 + t359 * t424 + (t262 + t427) * m(3) + t248;
t69 = (-t220 * t329 + (t111 * t320 + t105) * t235) * t408;
t68 = (-t218 * t332 + (t110 * t323 + t104) * t233) * t409;
t67 = (-t216 * t335 + (t109 * t326 + t103) * t231) * t410;
t66 = (-t108 * t220 + (-t105 * t406 - t114 * t363) * t235) * t139;
t65 = (-t107 * t218 + (-t104 * t406 - t113 * t364) * t233) * t138;
t64 = (-t106 * t216 + (-t103 * t406 - t112 * t365) * t231) * t137;
t60 = t63 ^ 2;
t59 = t62 ^ 2;
t58 = t61 ^ 2;
t57 = t136 * t63;
t56 = t135 * t62;
t55 = t134 * t61;
t54 = (-t177 * t265 + (-t222 * t93 + t81) * t235) * t408;
t53 = (-t176 * t266 + (-t222 * t91 + t79) * t233) * t409;
t52 = (-t175 * t267 + (-t222 * t89 + t77) * t231) * t410;
t51 = (t174 * t265 + (-t222 * t92 + t80) * t235) * t408;
t50 = (t173 * t266 + (-t222 * t90 + t78) * t233) * t409;
t49 = (t172 * t267 + (-t222 * t88 + t76) * t231) * t410;
t42 = (-t108 * t340 + (t114 * t93 - t81 * t406) * t235) * t139;
t41 = (-t107 * t343 + (t113 * t91 - t79 * t406) * t233) * t138;
t40 = (-t106 * t346 + (t112 * t89 - t77 * t406) * t231) * t137;
t39 = (t108 * t349 + (t114 * t92 - t80 * t406) * t235) * t139;
t38 = (t107 * t352 + (t113 * t90 - t78 * t406) * t233) * t138;
t37 = (t106 * t355 + (t112 * t88 - t76 * t406) * t231) * t137;
t33 = (-t84 * t340 + (t108 * t93 + t81 * t289) * t235) * t139;
t32 = (-t83 * t343 + (t107 * t91 + t79 * t290) * t233) * t138;
t31 = (-t82 * t346 + (t106 * t89 + t77 * t291) * t231) * t137;
t30 = (t84 * t349 + (t108 * t92 + t80 * t289) * t235) * t139;
t29 = (t83 * t352 + (t107 * t90 + t78 * t290) * t233) * t138;
t28 = (t82 * t355 + (t106 * t88 + t76 * t291) * t231) * t137;
t18 = (t27 * t432 + 0.2e1 * t283 + t366) * t102 * t139;
t17 = (t26 * t433 + 0.2e1 * t284 + t369) * t101 * t138;
t16 = (t25 * t434 + 0.2e1 * t285 + t372) * t100 * t137;
t15 = (-(t228 * t283 + t393 * t213 + (t160 * t219 * t431 + t149 + t270) * t102) * t219 * t102 + (-qJ(3,1) * t196 * t366 + (t27 + t277) * t311 + t27 * t160) * t63 + (t160 * t63 - t387) * t36) * t360;
t14 = (-(t228 * t284 + t394 * t211 + (t159 * t217 * t431 + t148 + t273) * t101) * t217 * t101 + (-qJ(3,2) * t195 * t369 + (t26 + t278) * t312 + t26 * t159) * t62 + (t159 * t62 - t388) * t35) * t361;
t13 = (-(t228 * t285 + t395 * t209 + (t158 * t215 * t431 + t147 + t276) * t100) * t215 * t100 + (-qJ(3,3) * t194 * t372 + (t25 + t279) * t313 + t25 * t158) * t61 + (t158 * t61 - t389) * t34) * t362;
t12 = ((-(t199 - 0.3e1 * t234) * t311 * t368 + (t228 * (-t431 * t63 + t36) * qJ(3,1) + (-0.3e1 * (-t234 / 0.3e1 + t199) * t171 + (t234 - t199) * t305) * t102) * t196 + (-t328 * t381 + ((-0.4e1 * t288 - t24) * t213 - t102 * (0.3e1 * t234 + t259)) * t229) * t219 - qJ(3,1) * (t149 * t367 + t393)) * t102 + ((t21 * t229 + t269 * t366) * t219 + t21 * pkin(1) + (t21 * t213 + (t417 + 0.1e1) * t102 * t314) * qJ(3,1)) * t63 + ((pkin(1) * t63 - t387) * t229 + (t63 * t327 + t412 * t366) * qJ(3,1)) * t36) * t360;
t11 = ((-(t199 - 0.3e1 * t232) * t312 * t371 + (t228 * (-t431 * t62 + t35) * qJ(3,2) + (-0.3e1 * (-t232 / 0.3e1 + t199) * t170 + (t232 - t199) * t305) * t101) * t195 + (-t331 * t382 + ((-0.4e1 * t287 - t23) * t211 - t101 * (0.3e1 * t232 + t259)) * t229) * t217 - qJ(3,2) * (t148 * t370 + t394)) * t101 + ((t20 * t229 + t272 * t369) * t217 + t20 * pkin(1) + (t20 * t211 + (t418 + 0.1e1) * t101 * t314) * qJ(3,2)) * t62 + ((pkin(1) * t62 - t388) * t229 + (t62 * t330 + t413 * t369) * qJ(3,2)) * t35) * t361;
t10 = ((-(t199 - 0.3e1 * t230) * t313 * t374 + (t228 * (-t431 * t61 + t34) * qJ(3,3) + (-0.3e1 * (-t230 / 0.3e1 + t199) * t169 + (t230 - t199) * t305) * t100) * t194 + (-t334 * t383 + ((-0.4e1 * t286 - t22) * t209 - t100 * (0.3e1 * t230 + t259)) * t229) * t215 - qJ(3,3) * (t147 * t373 + t395)) * t100 + ((t19 * t229 + t275 * t372) * t215 + t19 * pkin(1) + (t19 * t209 + (t419 + 0.1e1) * t100 * t314) * qJ(3,3)) * t61 + ((pkin(1) * t61 - t389) * t229 + (t61 * t333 + t414 * t372) * qJ(3,3)) * t34) * t362;
t9 = t15 * t406 - t18 * t289 + (-t12 + t60 * t202 + (-t412 * t202 + (-pkin(1) - t320) * t213) * t99 + t142 * t219 - t255 * t213) * m(3);
t8 = t14 * t406 - t17 * t290 + (-t11 + t59 * t201 + (-t413 * t201 + (-pkin(1) - t323) * t211) * t98 + t141 * t217 - t256 * t211) * m(3);
t7 = t13 * t406 - t16 * t291 + (-t10 + t58 * t200 + (-t414 * t200 + (-pkin(1) - t326) * t209) * t97 + t140 * t215 - t257 * t209) * m(3);
t6 = -t108 * t18 - t114 * t15 + (-t142 * t156 + t255 * t155) * t219 + (t142 * t155 + t255 * t156) * t213 + (t222 * t12 - 0.2e1 * t63 * t390) * m(3) + (t136 * t417 + (t117 * t213 + t400) * t219 + t156 * t397 + t136) * t99;
t5 = -t107 * t17 - t113 * t14 + (-t141 * t156 + t256 * t154) * t217 + (t141 * t154 + t256 * t156) * t211 + (t222 * t11 - 0.2e1 * t62 * t391) * m(3) + (t135 * t418 + (t116 * t211 + t401) * t217 + t156 * t398 + t135) * t98;
t4 = -t106 * t16 - t112 * t13 + (-t140 * t156 + t257 * t153) * t215 + (t140 * t153 + t257 * t156) * t209 + (t222 * t10 - 0.2e1 * t61 * t392) * m(3) + (t134 * t419 + (t115 * t209 + t402) * t215 + t156 * t399 + t134) * t97;
t3 = -t84 * t18 - t108 * t15 - t12 * t289 - 0.4e1 * (-t57 - t298 / 0.2e1) * t368 + (-0.2e1 * (t117 * t63 - t36 * t406) * t367 - t63 * (t130 * t63 + t36 * t304 + t400 * t420)) * t219 + ((t263 + t292) * t60 + (m(3) * t36 - t156 * t63) * t102 * t305) * t213 + (-t57 - t298) * t420 + (t251 * g(3) + t145 * t146 + t151) * t220 + (-t396 + (t152 + t251) * t145) * t214;
t2 = -t83 * t17 - t107 * t14 - t11 * t290 - 0.4e1 * (-t56 - t299 / 0.2e1) * t371 + (-0.2e1 * (t116 * t62 - t35 * t406) * t370 - t62 * (t130 * t62 + t35 * t304 + t401 * t421)) * t217 + ((t263 + t293) * t59 + (m(3) * t35 - t156 * t62) * t101 * t305) * t211 + (-t56 - t299) * t421 + (t252 * g(3) + t144 * t146 + t151) * t218 + (-t396 + (t152 + t252) * t144) * t212;
t1 = -t82 * t16 - t106 * t13 - t10 * t291 - 0.4e1 * (-t55 - t300 / 0.2e1) * t374 + (-0.2e1 * (t115 * t61 - t34 * t406) * t373 - t61 * (t130 * t61 + t34 * t304 + t402 * t422)) * t215 + ((t263 + t294) * t58 + (m(3) * t34 - t156 * t61) * t100 * t305) * t209 + (-t55 - t300) * t422 + (t253 * g(3) + t143 * t146 + t151) * t216 + (-t396 + (t152 + t253) * t143) * t210;
t43 = [(-g(1) + t208) * m(4) + (-t33 * t319 + (t33 * t350 + (-t208 * t33 - t3) * t177) * t214 + ((t42 * t93 + t54 * t81) * t208 + (t42 * t92 + t54 * t80) * t207 + (t105 * t54 - t42 * t363) * t206 + t93 * t6 + t81 * t9) * t235) * t139 + (-t32 * t322 + (t32 * t353 + (-t208 * t32 - t2) * t176) * t212 + ((t41 * t91 + t53 * t79) * t208 + (t41 * t90 + t53 * t78) * t207 + (t104 * t53 - t41 * t364) * t206 + t91 * t5 + t79 * t8) * t233) * t138 + (-t31 * t325 + (t31 * t356 + (-t208 * t31 - t1) * t175) * t210 + ((t40 * t89 + t52 * t77) * t208 + (t40 * t88 + t52 * t76) * t207 + (t103 * t52 - t40 * t365) * t206 + t89 * t4 + t77 * t7) * t231) * t137; (-g(2) + t207) * m(4) + (-t30 * t319 + (-t30 * t341 + (t207 * t30 + t3) * t174) * t214 + ((t39 * t93 + t51 * t81) * t208 + (t39 * t92 + t51 * t80) * t207 + (t105 * t51 - t39 * t363) * t206 + t92 * t6 + t80 * t9) * t235) * t139 + (-t29 * t322 + (-t29 * t344 + (t207 * t29 + t2) * t173) * t212 + ((t38 * t91 + t50 * t79) * t208 + (t38 * t90 + t50 * t78) * t207 + (t104 * t50 - t38 * t364) * t206 + t90 * t5 + t78 * t8) * t233) * t138 + (-t28 * t325 + (-t28 * t347 + (t207 * t28 + t1) * t172) * t210 + ((t37 * t89 + t49 * t77) * t208 + (t37 * t88 + t49 * t76) * t207 + (t103 * t49 - t37 * t365) * t206 + t88 * t4 + t76 * t7) * t231) * t137; (-g(3) + t206) * m(4) + (-t220 * t3 + (-t319 + (-t341 + t350) * t214) * (-t220 * t84 + (t105 * t289 - t108 * t363) * t235) * t139 + ((t66 * t93 + t69 * t81) * t208 + (t66 * t92 + t69 * t80) * t207 + (t105 * t69 - t66 * t363) * t206 - t6 * t363 + t105 * t9) * t235) * t139 + (-t218 * t2 + (-t322 + (-t344 + t353) * t212) * (-t218 * t83 + (t104 * t290 - t107 * t364) * t233) * t138 + ((t65 * t91 + t68 * t79) * t208 + (t65 * t90 + t68 * t78) * t207 + (t104 * t68 - t65 * t364) * t206 - t5 * t364 + t104 * t8) * t233) * t138 + (-t216 * t1 + (-t325 + (-t347 + t356) * t210) * (-t216 * t82 + (t103 * t291 - t106 * t365) * t231) * t137 + ((t64 * t89 + t67 * t77) * t208 + (t64 * t88 + t67 * t76) * t207 + (t103 * t67 - t64 * t365) * t206 - t4 * t365 + t103 * t7) * t231) * t137;];
tauX  = t43;
