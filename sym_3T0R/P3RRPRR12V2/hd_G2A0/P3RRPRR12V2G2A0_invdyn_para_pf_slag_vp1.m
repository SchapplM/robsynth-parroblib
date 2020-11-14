% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:19
% EndTime: 2020-08-06 19:20:31
% DurationCPUTime: 11.58s
% Computational Cost: add. (86913->640), mult. (102525->1046), div. (11898->6), fcn. (76671->18), ass. (0->382)
t210 = sin(qJ(2,3));
t432 = 0.2e1 * t210;
t212 = sin(qJ(2,2));
t431 = 0.2e1 * t212;
t214 = sin(qJ(2,1));
t430 = 0.2e1 * t214;
t162 = pkin(1) * t210 + qJ(3,3);
t211 = sin(qJ(1,3));
t229 = pkin(5) - pkin(6);
t165 = t211 * t229;
t170 = t210 * qJ(3,3);
t217 = cos(qJ(1,3));
t216 = cos(qJ(2,3));
t195 = t216 ^ 2;
t230 = pkin(2) + pkin(3);
t339 = (qJ(3,3) + t230) * (-qJ(3,3) + t230);
t280 = t195 * t339;
t314 = t230 * t216;
t103 = t217 * t280 + ((0.2e1 * t170 + pkin(1)) * t217 + t165) * t314 + qJ(3,3) * (t162 * t217 + t210 * t165);
t226 = xDP(3);
t227 = xDP(2);
t228 = xDP(1);
t159 = t170 + pkin(1);
t428 = t159 + t314;
t138 = 0.1e1 / t428;
t232 = 0.1e1 / qJ(3,3);
t363 = t138 * t232;
t264 = pkin(1) * t211 - t217 * t229;
t333 = t211 * qJ(3,3);
t121 = t264 * t210 + t333;
t284 = t210 * t333;
t124 = t264 + 0.2e1 * t284;
t279 = t210 * t339;
t401 = pkin(1) * qJ(3,3);
t131 = -t279 + t401;
t204 = legFrame(3,2);
t173 = sin(t204);
t176 = cos(t204);
t278 = t211 * t339;
t414 = -0.2e1 * t230;
t299 = qJ(3,3) * t414;
t346 = t176 * t230;
t355 = t173 * t230;
t377 = qJ(3,3) * t173;
t76 = (-t173 * t278 + t176 * t299) * t195 + (-t124 * t355 - t131 * t176) * t216 - t121 * t377 + t162 * t346;
t376 = qJ(3,3) * t176;
t77 = (t173 * t299 + t176 * t278) * t195 + (t124 * t346 - t131 * t173) * t216 + t121 * t376 + t162 * t355;
t34 = (t103 * t226 + t227 * t76 + t228 * t77) * t363;
t109 = t428 * t217 + t165;
t366 = t109 * t216;
t125 = t264 + t284;
t331 = t211 * t230;
t334 = t210 * t230;
t88 = (-t173 * t331 - t376) * t195 + (-t125 * t173 + t176 * t334) * t216 + t176 * t162;
t89 = (t176 * t331 - t377) * t195 + (t125 * t176 + t173 * t334) * t216 + t173 * t162;
t61 = (t226 * t366 + t227 * t88 + t228 * t89) * t363;
t387 = t230 * t61;
t25 = t34 - t387;
t163 = pkin(1) * t212 + qJ(3,2);
t213 = sin(qJ(1,2));
t166 = t213 * t229;
t171 = t212 * qJ(3,2);
t219 = cos(qJ(1,2));
t218 = cos(qJ(2,2));
t196 = t218 ^ 2;
t338 = (qJ(3,2) + t230) * (-qJ(3,2) + t230);
t277 = t196 * t338;
t313 = t230 * t218;
t104 = t219 * t277 + ((0.2e1 * t171 + pkin(1)) * t219 + t166) * t313 + qJ(3,2) * (t163 * t219 + t212 * t166);
t160 = t171 + pkin(1);
t427 = t160 + t313;
t139 = 0.1e1 / t427;
t234 = 0.1e1 / qJ(3,2);
t362 = t139 * t234;
t263 = pkin(1) * t213 - t219 * t229;
t327 = t213 * qJ(3,2);
t122 = t263 * t212 + t327;
t285 = t212 * t327;
t126 = t263 + 0.2e1 * t285;
t276 = t212 * t338;
t402 = pkin(1) * qJ(3,2);
t132 = -t276 + t402;
t205 = legFrame(2,2);
t174 = sin(t205);
t177 = cos(t205);
t275 = t213 * t338;
t300 = qJ(3,2) * t414;
t343 = t177 * t230;
t352 = t174 * t230;
t379 = qJ(3,2) * t174;
t78 = (-t174 * t275 + t177 * t300) * t196 + (-t126 * t352 - t132 * t177) * t218 - t122 * t379 + t163 * t343;
t378 = qJ(3,2) * t177;
t79 = (t174 * t300 + t177 * t275) * t196 + (t126 * t343 - t132 * t174) * t218 + t122 * t378 + t163 * t352;
t35 = (t104 * t226 + t227 * t78 + t228 * t79) * t362;
t110 = t427 * t219 + t166;
t365 = t110 * t218;
t127 = t263 + t285;
t325 = t213 * t230;
t328 = t212 * t230;
t90 = (-t174 * t325 - t378) * t196 + (-t127 * t174 + t177 * t328) * t218 + t177 * t163;
t91 = (t177 * t325 - t379) * t196 + (t127 * t177 + t174 * t328) * t218 + t174 * t163;
t62 = (t226 * t365 + t227 * t90 + t228 * t91) * t362;
t386 = t230 * t62;
t26 = t35 - t386;
t164 = pkin(1) * t214 + qJ(3,1);
t215 = sin(qJ(1,1));
t167 = t215 * t229;
t172 = t214 * qJ(3,1);
t221 = cos(qJ(1,1));
t220 = cos(qJ(2,1));
t197 = t220 ^ 2;
t337 = (qJ(3,1) + t230) * (-qJ(3,1) + t230);
t274 = t197 * t337;
t312 = t230 * t220;
t105 = t221 * t274 + ((0.2e1 * t172 + pkin(1)) * t221 + t167) * t312 + qJ(3,1) * (t164 * t221 + t214 * t167);
t161 = t172 + pkin(1);
t426 = t161 + t312;
t140 = 0.1e1 / t426;
t236 = 0.1e1 / qJ(3,1);
t361 = t140 * t236;
t262 = pkin(1) * t215 - t221 * t229;
t321 = t215 * qJ(3,1);
t123 = t262 * t214 + t321;
t286 = t214 * t321;
t128 = t262 + 0.2e1 * t286;
t273 = t214 * t337;
t403 = pkin(1) * qJ(3,1);
t133 = -t273 + t403;
t206 = legFrame(1,2);
t175 = sin(t206);
t178 = cos(t206);
t272 = t215 * t337;
t301 = qJ(3,1) * t414;
t340 = t178 * t230;
t349 = t175 * t230;
t381 = qJ(3,1) * t175;
t80 = (-t175 * t272 + t178 * t301) * t197 + (-t128 * t349 - t133 * t178) * t220 - t123 * t381 + t164 * t340;
t380 = qJ(3,1) * t178;
t81 = (t175 * t301 + t178 * t272) * t197 + (t128 * t340 - t133 * t175) * t220 + t123 * t380 + t164 * t349;
t36 = (t105 * t226 + t227 * t80 + t228 * t81) * t361;
t111 = t426 * t221 + t167;
t364 = t111 * t220;
t129 = t262 + t286;
t319 = t215 * t230;
t322 = t214 * t230;
t92 = (-t175 * t319 - t380) * t197 + (-t129 * t175 + t178 * t322) * t220 + t178 * t164;
t93 = (t178 * t319 - t381) * t197 + (t129 * t178 + t175 * t322) * t220 + t175 * t164;
t63 = (t226 * t364 + t227 * t92 + t228 * t93) * t361;
t385 = t230 * t63;
t27 = t36 - t385;
t429 = 0.2e1 * t230;
t231 = qJ(3,3) ^ 2;
t421 = 2 * rSges(3,3);
t425 = qJ(3,3) * t421 + t231;
t233 = qJ(3,2) ^ 2;
t424 = qJ(3,2) * t421 + t233;
t235 = qJ(3,1) ^ 2;
t423 = qJ(3,1) * t421 + t235;
t422 = -0.2e1 * pkin(1);
t100 = (-t211 * t226 + (-t173 * t227 + t176 * t228) * t217) * t138;
t420 = 0.2e1 * t100;
t101 = (-t213 * t226 + (-t174 * t227 + t177 * t228) * t219) * t139;
t419 = 0.2e1 * t101;
t102 = (-t215 * t226 + (-t175 * t227 + t178 * t228) * t221) * t140;
t418 = 0.2e1 * t102;
t417 = -0.2e1 * t195;
t416 = -0.2e1 * t196;
t415 = -0.2e1 * t197;
t413 = m(2) * rSges(2,1);
t225 = m(2) * rSges(2,2);
t224 = rSges(2,3) + pkin(5);
t223 = pkin(2) + rSges(3,1);
t412 = t195 - 0.1e1;
t411 = t196 - 0.1e1;
t410 = t197 - 0.1e1;
t409 = m(2) * t224;
t408 = m(3) * t138;
t407 = m(3) * t139;
t406 = m(3) * t140;
t222 = pkin(5) + rSges(3,2);
t405 = m(3) * t222;
t404 = m(3) * t223;
t201 = -qJ(3,3) - rSges(3,3);
t154 = m(3) * t201 + t225;
t400 = pkin(1) * t154;
t203 = -qJ(3,1) - rSges(3,3);
t156 = m(3) * t203 + t225;
t399 = pkin(1) * t156;
t202 = -qJ(3,2) - rSges(3,3);
t155 = m(3) * t202 + t225;
t398 = t155 * pkin(1);
t157 = t404 + t413;
t397 = t157 * pkin(1);
t22 = t229 * t25;
t290 = t100 * t401;
t396 = 0.2e1 * t290 + t22;
t23 = t229 * t26;
t291 = t101 * t402;
t395 = 0.2e1 * t291 + t23;
t24 = t229 * t27;
t292 = t102 * t403;
t394 = 0.2e1 * t292 + t24;
t393 = t201 * t34;
t392 = t202 * t35;
t391 = t203 * t36;
t335 = t210 * t229;
t283 = t100 * t335;
t390 = t216 * (t283 - t387);
t329 = t212 * t229;
t282 = t101 * t329;
t389 = t218 * (t282 - t386);
t323 = t214 * t229;
t281 = t102 * t323;
t388 = t220 * (t281 - t385);
t384 = t231 * t61;
t383 = t233 * t62;
t382 = t235 * t63;
t375 = t100 * t195;
t374 = t100 * t210;
t373 = t100 * t229;
t372 = t101 * t196;
t371 = t101 * t212;
t370 = t101 * t229;
t369 = t102 * t197;
t368 = t102 * t214;
t367 = t102 * t229;
t360 = t154 * t210;
t359 = t155 * t212;
t358 = t156 * t214;
t208 = xDDP(2);
t357 = t173 * t208;
t356 = t173 * t217;
t354 = t174 * t208;
t353 = t174 * t219;
t351 = t175 * t208;
t350 = t175 * t221;
t209 = xDDP(1);
t348 = t176 * t209;
t347 = t176 * t217;
t345 = t177 * t209;
t344 = t177 * t219;
t342 = t178 * t209;
t341 = t178 * t221;
t336 = t210 * t222;
t207 = xDDP(3);
t332 = t211 * t207;
t330 = t212 * t222;
t326 = t213 * t207;
t324 = t214 * t222;
t320 = t215 * t207;
t318 = t216 * t223;
t317 = t218 * t223;
t316 = t220 * t223;
t315 = t229 * t230;
t311 = -Icges(2,1) - Icges(3,1);
t310 = Icges(2,6) - Icges(3,6);
t309 = 0.2e1 * pkin(1);
t308 = -0.2e1 * t405;
t307 = rSges(3,3) - t223;
t306 = rSges(3,3) + t223;
t304 = m(3) * t393;
t303 = m(3) * t392;
t302 = m(3) * t391;
t298 = t201 * t405;
t297 = t202 * t405;
t296 = t203 * t405;
t295 = m(3) * t336;
t294 = m(3) * t330;
t293 = m(3) * t324;
t289 = t216 * qJ(3,3) * t61;
t288 = t218 * qJ(3,2) * t62;
t287 = t63 * t220 * qJ(3,1);
t271 = t217 * t336;
t270 = t219 * t330;
t269 = t221 * t324;
t268 = -rSges(2,1) * t225 + Icges(2,4) - Icges(3,5);
t267 = rSges(2,2) * t409 - t310;
t237 = rSges(3,3) ^ 2;
t244 = pkin(1) ^ 2;
t266 = t222 ^ 2 + t237 + t244;
t261 = pkin(2) ^ 2 + t237 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t260 = (pkin(6) ^ 2) + t244 + ((-2 * pkin(6) + pkin(5)) * pkin(5));
t238 = rSges(2,2) ^ 2;
t240 = rSges(2,1) ^ 2;
t259 = (t238 + t240) * m(2) + Icges(3,2) + Icges(2,3);
t144 = g(1) * t176 - g(2) * t173;
t258 = g(3) * t217 + t144 * t211;
t145 = g(1) * t177 - g(2) * t174;
t257 = g(3) * t219 + t145 * t213;
t146 = g(1) * t178 - g(2) * t175;
t256 = g(3) * t221 + t146 * t215;
t19 = t230 * t25 - t384;
t20 = t230 * t26 - t383;
t21 = t230 * t27 - t382;
t255 = Icges(2,2) + Icges(3,3) + (-t238 + t240) * m(2) + t311;
t254 = t157 * t216 - t360;
t253 = t157 * t218 - t359;
t252 = t157 * t220 - t358;
t251 = t222 * t404 - Icges(3,4) - Icges(2,5);
t250 = -t224 * t225 + t310;
t249 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (t224 ^ 2 + t238 + t244) * m(2) + Icges(1,3) - t311;
t248 = t224 * t413 + t251;
t200 = t230 ^ 2;
t153 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t152 = t153 * g(3);
t151 = 0.2e1 * t397;
t150 = t235 + t260;
t149 = t233 + t260;
t148 = t231 + t260;
t147 = m(1) * rSges(1,2) - t405 - t409;
t143 = g(1) * t175 + g(2) * t178;
t142 = g(1) * t174 + g(2) * t177;
t141 = g(1) * t173 + g(2) * t176;
t137 = t147 * g(3);
t136 = -t203 * t404 + t268;
t135 = -t202 * t404 + t268;
t134 = -t201 * t404 + t268;
t130 = rSges(2,1) * t409 + t251;
t117 = -(qJ(3,1) + t306) * (qJ(3,1) + t307) * m(3) + t255;
t116 = -(qJ(3,2) + t306) * (qJ(3,2) + t307) * m(3) + t255;
t115 = -(qJ(3,3) + t306) * (qJ(3,3) + t307) * m(3) + t255;
t114 = (t261 + t423) * m(3) + t259;
t113 = (t261 + t424) * m(3) + t259;
t112 = (t261 + t425) * m(3) + t259;
t108 = (t250 - t296) * t220 - t214 * t248;
t107 = (t250 - t297) * t218 - t212 * t248;
t106 = (t250 - t298) * t216 - t210 * t248;
t99 = t102 ^ 2;
t98 = t101 ^ 2;
t97 = t100 ^ 2;
t84 = t117 * t197 + (t136 * t430 + t151) * t220 + t358 * t422 + (t266 + t423) * m(3) + t249;
t83 = t116 * t196 + (t135 * t431 + t151) * t218 + t359 * t422 + (t266 + t424) * m(3) + t249;
t82 = t115 * t195 + (t134 * t432 + t151) * t216 + t360 * t422 + (t266 + t425) * m(3) + t249;
t69 = (-t215 * t324 + (-t111 * t316 + t105) * t236) * t406;
t68 = (-t213 * t330 + (-t110 * t317 + t104) * t234) * t407;
t67 = (-t211 * t336 + (-t109 * t318 + t103) * t232) * t408;
t66 = (-t108 * t215 + (-t105 * t404 + t114 * t364) * t236) * t140;
t65 = (-t107 * t213 + (-t104 * t404 + t113 * t365) * t234) * t139;
t64 = (-t106 * t211 + (-t103 * t404 + t112 * t366) * t232) * t138;
t60 = t63 ^ 2;
t59 = t62 ^ 2;
t58 = t61 ^ 2;
t57 = t136 * t63;
t56 = t135 * t62;
t55 = t134 * t61;
t54 = (t178 * t269 + (-t223 * t93 + t81) * t236) * t406;
t53 = (t177 * t270 + (-t223 * t91 + t79) * t234) * t407;
t52 = (t176 * t271 + (-t223 * t89 + t77) * t232) * t408;
t51 = (-t175 * t269 + (-t223 * t92 + t80) * t236) * t406;
t50 = (-t174 * t270 + (-t223 * t90 + t78) * t234) * t407;
t49 = (-t173 * t271 + (-t223 * t88 + t76) * t232) * t408;
t42 = (t108 * t341 + (t114 * t93 - t81 * t404) * t236) * t140;
t41 = (t107 * t344 + (t113 * t91 - t79 * t404) * t234) * t139;
t40 = (t106 * t347 + (t112 * t89 - t77 * t404) * t232) * t138;
t39 = (-t108 * t350 + (t114 * t92 - t80 * t404) * t236) * t140;
t38 = (-t107 * t353 + (t113 * t90 - t78 * t404) * t234) * t139;
t37 = (-t106 * t356 + (t112 * t88 - t76 * t404) * t232) * t138;
t33 = (t84 * t341 + (t108 * t93 + t81 * t293) * t236) * t140;
t32 = (t83 * t344 + (t107 * t91 + t79 * t294) * t234) * t139;
t31 = (t82 * t347 + (t106 * t89 + t77 * t295) * t232) * t138;
t30 = (-t84 * t350 + (t108 * t92 + t80 * t293) * t236) * t140;
t29 = (-t83 * t353 + (t107 * t90 + t78 * t294) * t234) * t139;
t28 = (-t82 * t356 + (t106 * t88 + t76 * t295) * t232) * t138;
t18 = (t27 * t430 + 0.2e1 * t287 + t367) * t102 * t140;
t17 = (t26 * t431 + 0.2e1 * t288 + t370) * t101 * t139;
t16 = (t25 * t432 + 0.2e1 * t289 + t373) * t100 * t138;
t15 = (-(t229 * t287 + t394 * t214 + (t161 * t220 * t429 + t150 + t274) * t102) * t220 * t102 + (-qJ(3,1) * t197 * t367 + (t27 + t281) * t312 + t27 * t161) * t63 + (t161 * t63 - t388) * t36) * t361;
t14 = (-(t229 * t288 + t395 * t212 + (t160 * t218 * t429 + t149 + t277) * t101) * t218 * t101 + (-qJ(3,2) * t196 * t370 + (t26 + t282) * t313 + t26 * t160) * t62 + (t160 * t62 - t389) * t35) * t362;
t13 = (-(t229 * t289 + t396 * t210 + (t159 * t216 * t429 + t148 + t280) * t100) * t216 * t100 + (-qJ(3,3) * t195 * t373 + (t25 + t283) * t314 + t25 * t159) * t61 + (t159 * t61 - t390) * t34) * t363;
t12 = ((-(t200 - 0.3e1 * t235) * t312 * t369 + (t229 * (-t429 * t63 + t36) * qJ(3,1) + (-0.3e1 * (-t235 / 0.3e1 + t200) * t172 + (t235 - t200) * t309) * t102) * t197 + (-t323 * t382 + ((-0.4e1 * t292 - t24) * t214 - t102 * (0.3e1 * t235 + t260)) * t230) * t220 - (t150 * t368 + t394) * qJ(3,1)) * t102 + ((t21 * t230 + t273 * t367) * t220 + t21 * pkin(1) + (t21 * t214 + (t415 + 0.1e1) * t102 * t315) * qJ(3,1)) * t63 + ((pkin(1) * t63 - t388) * t230 + (t63 * t322 + t410 * t367) * qJ(3,1)) * t36) * t361;
t11 = ((-(t200 - 0.3e1 * t233) * t313 * t372 + (t229 * (-t429 * t62 + t35) * qJ(3,2) + (-0.3e1 * (-t233 / 0.3e1 + t200) * t171 + (t233 - t200) * t309) * t101) * t196 + (-t329 * t383 + ((-0.4e1 * t291 - t23) * t212 - t101 * (0.3e1 * t233 + t260)) * t230) * t218 - (t149 * t371 + t395) * qJ(3,2)) * t101 + ((t20 * t230 + t276 * t370) * t218 + t20 * pkin(1) + (t20 * t212 + (t416 + 0.1e1) * t101 * t315) * qJ(3,2)) * t62 + ((pkin(1) * t62 - t389) * t230 + (t62 * t328 + t411 * t370) * qJ(3,2)) * t35) * t362;
t10 = ((-(t200 - 0.3e1 * t231) * t314 * t375 + (t229 * (-t429 * t61 + t34) * qJ(3,3) + (-0.3e1 * (-t231 / 0.3e1 + t200) * t170 + (t231 - t200) * t309) * t100) * t195 + (-t335 * t384 + ((-0.4e1 * t290 - t22) * t210 - t100 * (0.3e1 * t231 + t260)) * t230) * t216 - (t148 * t374 + t396) * qJ(3,3)) * t100 + ((t19 * t230 + t279 * t373) * t216 + t19 * pkin(1) + (t19 * t210 + (t417 + 0.1e1) * t100 * t315) * qJ(3,3)) * t61 + ((pkin(1) * t61 - t390) * t230 + (t61 * t334 + t412 * t373) * qJ(3,3)) * t34) * t363;
t9 = t15 * t404 - t18 * t293 + (-t12 + t60 * t203 + (-t410 * t203 + (-pkin(1) - t316) * t214) * t99 + t143 * t220 - t256 * t214) * m(3);
t8 = t14 * t404 - t17 * t294 + (-t11 + t59 * t202 + (-t411 * t202 + (-pkin(1) - t317) * t212) * t98 + t142 * t218 - t257 * t212) * m(3);
t7 = t13 * t404 - t16 * t295 + (-t10 + t58 * t201 + (-t412 * t201 + (-pkin(1) - t318) * t210) * t97 + t141 * t216 - t258 * t210) * m(3);
t6 = -t108 * t18 - t114 * t15 + (-t143 * t157 + t256 * t156) * t220 + (t143 * t156 + t256 * t157) * t214 + (t223 * t12 - 0.2e1 * t63 * t391) * m(3) + (t136 * t415 + (t117 * t214 + t399) * t220 + t214 * t397 + t136) * t99;
t5 = -t107 * t17 - t113 * t14 + (-t142 * t157 + t257 * t155) * t218 + (t142 * t155 + t257 * t157) * t212 + (t223 * t11 - 0.2e1 * t62 * t392) * m(3) + (t135 * t416 + (t116 * t212 + t398) * t218 + t212 * t397 + t135) * t98;
t4 = -t106 * t16 - t112 * t13 + (-t141 * t157 + t258 * t154) * t216 + (t141 * t154 + t258 * t157) * t210 + (t223 * t10 - 0.2e1 * t61 * t393) * m(3) + (t134 * t417 + (t115 * t210 + t400) * t216 + t210 * t397 + t134) * t97;
t3 = -t84 * t18 - t108 * t15 - t12 * t293 - 0.4e1 * (-t57 - t302 / 0.2e1) * t369 + (-0.2e1 * (t117 * t63 - t36 * t404) * t368 - t63 * (t130 * t63 + t36 * t308 + t399 * t418)) * t220 + ((t267 + t296) * t60 + (m(3) * t36 - t157 * t63) * t102 * t309) * t214 + (-t57 - t302) * t418 + (t137 + (-t153 - t252) * t146) * t221 + (t252 * g(3) + t146 * t147 + t152) * t215;
t2 = -t83 * t17 - t107 * t14 - t11 * t294 - 0.4e1 * (-t56 - t303 / 0.2e1) * t372 + (-0.2e1 * (t116 * t62 - t35 * t404) * t371 - t62 * (t130 * t62 + t35 * t308 + t398 * t419)) * t218 + ((t267 + t297) * t59 + (m(3) * t35 - t157 * t62) * t101 * t309) * t212 + (-t56 - t303) * t419 + (t137 + (-t153 - t253) * t145) * t219 + (t253 * g(3) + t145 * t147 + t152) * t213;
t1 = -t82 * t16 - t106 * t13 - t10 * t295 - 0.4e1 * (-t55 - t304 / 0.2e1) * t375 + (-0.2e1 * (t115 * t61 - t34 * t404) * t374 - t61 * (t130 * t61 + t34 * t308 + t400 * t420)) * t216 + ((t267 + t298) * t58 + (m(3) * t34 - t157 * t61) * t100 * t309) * t210 + (-t55 - t304) * t420 + (t137 + (-t153 - t254) * t144) * t217 + (t254 * g(3) + t144 * t147 + t152) * t211;
t43 = [(-g(1) + t209) * m(4) + (-t33 * t320 + (-t33 * t351 + (t209 * t33 + t3) * t178) * t221 + ((t42 * t93 + t54 * t81) * t209 + (t42 * t92 + t54 * t80) * t208 + (t105 * t54 + t42 * t364) * t207 + t93 * t6 + t81 * t9) * t236) * t140 + (-t32 * t326 + (-t32 * t354 + (t209 * t32 + t2) * t177) * t219 + ((t41 * t91 + t53 * t79) * t209 + (t41 * t90 + t53 * t78) * t208 + (t104 * t53 + t41 * t365) * t207 + t91 * t5 + t79 * t8) * t234) * t139 + (-t31 * t332 + (-t31 * t357 + (t209 * t31 + t1) * t176) * t217 + ((t40 * t89 + t52 * t77) * t209 + (t40 * t88 + t52 * t76) * t208 + (t103 * t52 + t40 * t366) * t207 + t89 * t4 + t77 * t7) * t232) * t138; (-g(2) + t208) * m(4) + (-t30 * t320 + (t30 * t342 + (-t208 * t30 - t3) * t175) * t221 + ((t39 * t93 + t51 * t81) * t209 + (t39 * t92 + t51 * t80) * t208 + (t105 * t51 + t39 * t364) * t207 + t92 * t6 + t80 * t9) * t236) * t140 + (-t29 * t326 + (t29 * t345 + (-t208 * t29 - t2) * t174) * t219 + ((t38 * t91 + t50 * t79) * t209 + (t38 * t90 + t50 * t78) * t208 + (t104 * t50 + t38 * t365) * t207 + t90 * t5 + t78 * t8) * t234) * t139 + (-t28 * t332 + (t28 * t348 + (-t208 * t28 - t1) * t173) * t217 + ((t37 * t89 + t49 * t77) * t209 + (t37 * t88 + t49 * t76) * t208 + (t103 * t49 + t37 * t366) * t207 + t88 * t4 + t76 * t7) * t232) * t138; (-g(3) + t207) * m(4) + (-t215 * t3 + (-t320 + (t342 - t351) * t221) * (-t215 * t84 + (t105 * t293 + t108 * t364) * t236) * t140 + ((t66 * t93 + t69 * t81) * t209 + (t66 * t92 + t69 * t80) * t208 + (t105 * t69 + t66 * t364) * t207 + t6 * t364 + t105 * t9) * t236) * t140 + (-t213 * t2 + (-t326 + (t345 - t354) * t219) * (-t213 * t83 + (t104 * t294 + t107 * t365) * t234) * t139 + ((t65 * t91 + t68 * t79) * t209 + (t65 * t90 + t68 * t78) * t208 + (t104 * t68 + t65 * t365) * t207 + t5 * t365 + t104 * t8) * t234) * t139 + (-t211 * t1 + (-t332 + (t348 - t357) * t217) * (-t211 * t82 + (t103 * t295 + t106 * t366) * t232) * t138 + ((t64 * t89 + t67 * t77) * t209 + (t64 * t88 + t67 * t76) * t208 + (t103 * t67 + t64 * t366) * t207 + t4 * t366 + t103 * t7) * t232) * t138;];
tauX  = t43;
