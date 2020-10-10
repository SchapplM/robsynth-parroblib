% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V2G1A0
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:14:29
% EndTime: 2020-08-06 19:14:40
% DurationCPUTime: 11.29s
% Computational Cost: add. (65157->615), mult. (78588->997), div. (9303->6), fcn. (67950->18), ass. (0->372)
t234 = pkin(5) - pkin(6);
t235 = pkin(2) + pkin(3);
t427 = t234 * t235;
t317 = 2 * pkin(1);
t426 = -0.2e1 * t235 * (-pkin(5) / 0.2e1 + pkin(6) / 0.2e1);
t221 = cos(qJ(2,3));
t198 = t221 ^ 2;
t209 = legFrame(3,3);
t178 = sin(t209);
t181 = cos(t209);
t216 = sin(qJ(1,3));
t222 = cos(qJ(1,3));
t118 = t178 * t216 - t181 * t222;
t119 = t178 * t222 + t181 * t216;
t215 = sin(qJ(2,3));
t354 = qJ(3,3) * t215;
t162 = pkin(1) + t354;
t323 = t235 * t221;
t143 = 0.1e1 / (t162 + t323);
t232 = xDP(2);
t233 = xDP(1);
t94 = (-t118 * t232 - t119 * t233) * t143;
t425 = t198 * t94;
t223 = cos(qJ(2,2));
t200 = t223 ^ 2;
t210 = legFrame(2,3);
t179 = sin(t210);
t182 = cos(t210);
t218 = sin(qJ(1,2));
t224 = cos(qJ(1,2));
t120 = t179 * t218 - t182 * t224;
t121 = t179 * t224 + t182 * t218;
t217 = sin(qJ(2,2));
t356 = qJ(3,2) * t217;
t164 = pkin(1) + t356;
t322 = t235 * t223;
t144 = 0.1e1 / (t164 + t322);
t95 = (-t120 * t232 - t121 * t233) * t144;
t424 = t200 * t95;
t225 = cos(qJ(2,1));
t202 = t225 ^ 2;
t211 = legFrame(1,3);
t180 = sin(t211);
t183 = cos(t211);
t220 = sin(qJ(1,1));
t226 = cos(qJ(1,1));
t122 = t180 * t220 - t183 * t226;
t123 = t180 * t226 + t183 * t220;
t219 = sin(qJ(2,1));
t358 = qJ(3,1) * t219;
t166 = pkin(1) + t358;
t321 = t235 * t225;
t145 = 0.1e1 / (t166 + t321);
t96 = (-t122 * t232 - t123 * t233) * t145;
t423 = t202 * t96;
t231 = xDP(3);
t237 = 0.1e1 / qJ(3,3);
t174 = t234 * t222;
t125 = t162 * t216 - t174;
t170 = t216 * t234;
t275 = t162 * t222 + t170;
t85 = t118 * t323 + t125 * t178 - t275 * t181;
t86 = t119 * t323 + t125 * t181 + t178 * t275;
t64 = (t215 * t231 + (t232 * t86 - t233 * t85) * t221 * t143) * t237;
t34 = t234 * t64;
t422 = t235 * t34;
t239 = 0.1e1 / qJ(3,2);
t175 = t234 * t224;
t127 = t164 * t218 - t175;
t171 = t218 * t234;
t273 = t164 * t224 + t171;
t87 = t120 * t322 + t127 * t179 - t273 * t182;
t88 = t121 * t322 + t127 * t182 + t179 * t273;
t65 = (t217 * t231 + (t232 * t88 - t233 * t87) * t223 * t144) * t239;
t35 = t234 * t65;
t421 = t235 * t35;
t241 = 0.1e1 / qJ(3,1);
t176 = t234 * t226;
t129 = t166 * t220 - t176;
t172 = t220 * t234;
t271 = t166 * t226 + t172;
t89 = t122 * t321 + t129 * t180 - t271 * t183;
t90 = t123 * t321 + t129 * t183 + t180 * t271;
t66 = (t219 * t231 + (t232 * t90 - t233 * t89) * t225 * t145) * t241;
t36 = t234 * t66;
t420 = t235 * t36;
t367 = t235 * t64;
t366 = t235 * t65;
t365 = t235 * t66;
t340 = t143 * t237;
t339 = t144 * t239;
t338 = t145 * t241;
t236 = qJ(3,3) ^ 2;
t415 = 2 * rSges(3,3);
t419 = qJ(3,3) * t415 + t236;
t238 = qJ(3,2) ^ 2;
t418 = qJ(3,2) * t415 + t238;
t240 = qJ(3,1) ^ 2;
t417 = qJ(3,1) * t415 + t240;
t414 = 0.2e1 * t94;
t413 = 0.2e1 * t95;
t412 = 0.2e1 * t96;
t411 = -0.2e1 * t198;
t410 = -0.2e1 * t200;
t409 = -0.2e1 * t202;
t408 = -0.3e1 * t236;
t407 = -0.3e1 * t238;
t406 = -0.3e1 * t240;
t405 = m(2) * rSges(2,1);
t230 = m(2) * rSges(2,2);
t229 = rSges(2,3) + pkin(5);
t228 = pkin(2) + rSges(3,1);
t404 = t198 - 0.1e1;
t403 = t200 - 0.1e1;
t402 = t202 - 0.1e1;
t401 = m(2) * t229;
t400 = m(3) * t143;
t399 = m(3) * t144;
t398 = m(3) * t145;
t227 = pkin(5) + rSges(3,2);
t397 = m(3) * t227;
t396 = m(3) * t228;
t206 = -qJ(3,3) - rSges(3,3);
t153 = m(3) * t206 + t230;
t395 = pkin(1) * t153;
t207 = -qJ(3,2) - rSges(3,3);
t154 = m(3) * t207 + t230;
t394 = pkin(1) * t154;
t208 = -qJ(3,1) - rSges(3,3);
t155 = m(3) * t208 + t230;
t393 = pkin(1) * t155;
t159 = t396 + t405;
t392 = g(3) * t159;
t391 = t215 * pkin(1);
t390 = t217 * pkin(1);
t389 = t219 * pkin(1);
t326 = t235 * t215;
t353 = qJ(3,3) * t221;
t147 = t326 - t353;
t167 = qJ(3,3) + t391;
t115 = t167 * t216 - t215 * t174;
t161 = pkin(1) + 0.2e1 * t354;
t124 = t161 * t216 - t174;
t262 = t167 * t222 + t215 * t170;
t276 = t161 * t222 + t170;
t336 = (qJ(3,3) + t235) * (-qJ(3,3) + t235);
t282 = t198 * t336;
t73 = -t118 * t282 - (t178 * t124 - t276 * t181) * t323 - (t178 * t115 - t262 * t181) * qJ(3,3);
t74 = t119 * t282 + (t124 * t181 + t276 * t178) * t323 + (t115 * t181 + t262 * t178) * qJ(3,3);
t46 = t147 * t237 * t231 + (t232 * t74 + t233 * t73) * t340;
t388 = t206 * t46;
t325 = t235 * t217;
t355 = qJ(3,2) * t223;
t148 = t325 - t355;
t168 = qJ(3,2) + t390;
t116 = t168 * t218 - t217 * t175;
t163 = pkin(1) + 0.2e1 * t356;
t126 = t163 * t218 - t175;
t261 = t168 * t224 + t217 * t171;
t274 = t163 * t224 + t171;
t335 = (qJ(3,2) + t235) * (-qJ(3,2) + t235);
t281 = t200 * t335;
t75 = -t120 * t281 - (t179 * t126 - t274 * t182) * t322 - (t179 * t116 - t261 * t182) * qJ(3,2);
t76 = t121 * t281 + (t126 * t182 + t274 * t179) * t322 + (t116 * t182 + t261 * t179) * qJ(3,2);
t47 = t148 * t239 * t231 + (t232 * t76 + t233 * t75) * t339;
t387 = t207 * t47;
t324 = t235 * t219;
t357 = qJ(3,1) * t225;
t149 = t324 - t357;
t169 = qJ(3,1) + t389;
t117 = t169 * t220 - t219 * t176;
t165 = pkin(1) + 0.2e1 * t358;
t128 = t165 * t220 - t176;
t260 = t169 * t226 + t219 * t172;
t272 = t165 * t226 + t172;
t334 = (qJ(3,1) + t235) * (-qJ(3,1) + t235);
t280 = t202 * t334;
t77 = -t122 * t280 - (t180 * t128 - t272 * t183) * t321 - (t180 * t117 - t260 * t183) * qJ(3,1);
t78 = t123 * t280 + (t128 * t183 + t272 * t180) * t321 + (t117 * t183 + t260 * t180) * qJ(3,1);
t48 = t149 * t241 * t231 + (t232 * t78 + t233 * t77) * t338;
t386 = t208 * t48;
t385 = t215 * t94;
t384 = t217 * t95;
t383 = t219 * t96;
t361 = t94 * t234;
t288 = t215 * t361;
t382 = t221 * (t288 - t367);
t381 = t221 * t85;
t380 = t221 * t86;
t243 = rSges(2,2) ^ 2;
t245 = rSges(2,1) ^ 2;
t269 = (t243 + t245) * m(2) + Icges(3,2) + Icges(2,3);
t242 = rSges(3,3) ^ 2;
t270 = pkin(2) ^ 2 + t242 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t103 = (t270 + t419) * m(3) + t269;
t97 = (t103 * t215 - t147 * t396) * t237;
t379 = t221 * t97;
t360 = t95 * t234;
t287 = t217 * t360;
t378 = t223 * (t287 - t366);
t377 = t223 * t87;
t376 = t223 * t88;
t104 = (t270 + t418) * m(3) + t269;
t98 = (t104 * t217 - t148 * t396) * t239;
t375 = t223 * t98;
t359 = t96 * t234;
t286 = t219 * t359;
t374 = t225 * (t286 - t365);
t373 = t225 * t89;
t372 = t225 * t90;
t105 = (t270 + t417) * m(3) + t269;
t99 = (t105 * t219 - t149 * t396) * t241;
t371 = t225 * t99;
t364 = t235 * t94;
t363 = t235 * t95;
t362 = t235 * t96;
t263 = t227 * t396 - Icges(3,4) - Icges(2,5);
t254 = t229 * t405 + t263;
t319 = Icges(2,6) - Icges(3,6);
t256 = -t229 * t230 + t319;
t303 = t206 * t397;
t100 = (t256 - t303) * t221 - t215 * t254;
t352 = t100 * t221;
t302 = t207 * t397;
t101 = (t256 - t302) * t223 - t217 * t254;
t351 = t101 * t223;
t301 = t208 * t397;
t102 = (t256 - t301) * t225 - t219 * t254;
t350 = t102 * t225;
t349 = t103 * t221;
t348 = t104 * t223;
t347 = t105 * t225;
t213 = xDDP(2);
t346 = t118 * t213;
t214 = xDDP(1);
t345 = t119 * t214;
t344 = t120 * t213;
t343 = t121 * t214;
t342 = t122 * t213;
t341 = t123 * t214;
t247 = pkin(1) ^ 2;
t152 = pkin(6) ^ 2 + t247 + (-0.2e1 * pkin(6) + pkin(5)) * pkin(5);
t337 = t152 * t235;
t333 = t215 * t227;
t332 = t217 * t227;
t331 = t219 * t227;
t330 = t221 * t228;
t329 = t223 * t228;
t328 = t225 * t228;
t320 = -Icges(2,1) - Icges(3,1);
t318 = -t234 ^ 2 - t247;
t316 = -0.2e1 * t397;
t315 = -0.2e1 * t391;
t314 = -0.2e1 * t390;
t313 = -0.2e1 * t389;
t312 = rSges(3,3) - t228;
t311 = rSges(3,3) + t228;
t310 = 0.2e1 * t235;
t309 = m(3) * t388;
t308 = m(3) * t387;
t307 = m(3) * t386;
t306 = pkin(1) * qJ(3,1) * t96;
t305 = pkin(1) * qJ(3,2) * t95;
t304 = pkin(1) * qJ(3,3) * t94;
t300 = m(3) * t333;
t299 = m(3) * t332;
t298 = m(3) * t331;
t297 = t46 * t340;
t296 = t64 * t340;
t295 = t94 * t340;
t294 = t47 * t339;
t293 = t65 * t339;
t292 = t95 * t339;
t291 = t48 * t338;
t290 = t66 * t338;
t289 = t96 * t338;
t285 = qJ(3,1) * t202 * t234;
t284 = qJ(3,2) * t200 * t234;
t283 = qJ(3,3) * t198 * t234;
t279 = -rSges(2,1) * t230 + Icges(2,4) - Icges(3,5);
t278 = rSges(2,2) * t401 - t319;
t277 = t227 ^ 2 + t242 + t247;
t22 = t47 - t366;
t20 = t235 * t22 - t238 * t65;
t23 = t48 - t365;
t21 = t235 * t23 - t240 * t66;
t24 = t46 - t367;
t19 = t235 * t24 - t236 * t64;
t268 = Icges(2,2) + Icges(3,3) + (-t243 + t245) * m(2) + t320;
t137 = -g(1) * t178 + g(2) * t181;
t140 = g(1) * t181 + g(2) * t178;
t266 = t137 * t216 + t140 * t222;
t138 = -g(1) * t179 + g(2) * t182;
t141 = g(1) * t182 + g(2) * t179;
t265 = t138 * t218 + t141 * t224;
t139 = -g(1) * t180 + g(2) * t183;
t142 = g(1) * t183 + g(2) * t180;
t264 = t139 * t220 + t142 * t226;
t151 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t259 = t153 * t215 - t159 * t221 - t151;
t258 = t154 * t217 - t159 * t223 - t151;
t257 = t155 * t219 - t159 * t225 - t151;
t255 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (t229 ^ 2 + t243 + t247) * m(2) + Icges(1,3) - t320;
t212 = xDDP(3);
t205 = t235 ^ 2;
t201 = t225 * t202;
t199 = t223 * t200;
t197 = t221 * t198;
t150 = t159 * t317;
t146 = m(1) * rSges(1,2) - t397 - t401;
t136 = -t208 * t396 + t279;
t135 = -t207 * t396 + t279;
t134 = -t206 * t396 + t279;
t133 = rSges(2,1) * t401 + t263;
t111 = (-t219 * t228 + t149) * t241 * m(3);
t110 = (-t217 * t228 + t148) * t239 * m(3);
t109 = (-t215 * t228 + t147) * t237 * m(3);
t108 = -(qJ(3,1) + t311) * (qJ(3,1) + t312) * m(3) + t268;
t107 = -(qJ(3,2) + t311) * (qJ(3,2) + t312) * m(3) + t268;
t106 = -(qJ(3,3) + t311) * (qJ(3,3) + t312) * m(3) + t268;
t93 = t96 ^ 2;
t92 = t95 ^ 2;
t91 = t94 ^ 2;
t81 = t108 * t202 + (0.2e1 * t136 * t219 + t150) * t225 + t155 * t313 + (t277 + t417) * m(3) + t255;
t80 = t107 * t200 + (0.2e1 * t135 * t217 + t150) * t223 + t154 * t314 + (t277 + t418) * m(3) + t255;
t79 = t106 * t198 + (0.2e1 * t134 * t215 + t150) * t221 + t153 * t315 + (t277 + t419) * m(3) + t255;
t63 = t66 ^ 2;
t62 = t65 ^ 2;
t61 = t64 ^ 2;
t60 = t136 * t66;
t59 = t135 * t65;
t58 = t134 * t64;
t54 = (-t123 * t331 + (t89 * t328 + t77) * t241) * t398;
t53 = (-t122 * t331 + (-t90 * t328 + t78) * t241) * t398;
t52 = (-t121 * t332 + (t87 * t329 + t75) * t239) * t399;
t51 = (-t120 * t332 + (-t88 * t329 + t76) * t239) * t399;
t50 = (-t119 * t333 + (t85 * t330 + t73) * t237) * t400;
t49 = (-t118 * t333 + (-t86 * t330 + t74) * t237) * t400;
t45 = t48 * t234;
t44 = t47 * t234;
t43 = t46 * t234;
t42 = (-t102 * t123 + (-t89 * t347 - t77 * t396) * t241) * t145;
t41 = (-t102 * t122 + (t90 * t347 - t78 * t396) * t241) * t145;
t40 = (-t101 * t121 + (-t87 * t348 - t75 * t396) * t239) * t144;
t39 = (-t101 * t120 + (t88 * t348 - t76 * t396) * t239) * t144;
t38 = (-t100 * t119 + (-t85 * t349 - t73 * t396) * t237) * t143;
t37 = (-t100 * t118 + (t86 * t349 - t74 * t396) * t237) * t143;
t30 = (-t123 * t81 + (t77 * t298 - t89 * t350) * t241) * t145;
t29 = (-t122 * t81 + (t78 * t298 + t90 * t350) * t241) * t145;
t28 = (-t121 * t80 + (t75 * t299 - t87 * t351) * t239) * t144;
t27 = (-t120 * t80 + (t76 * t299 + t88 * t351) * t239) * t144;
t26 = (-t119 * t79 + (t73 * t300 - t85 * t352) * t237) * t143;
t25 = (-t118 * t79 + (t74 * t300 + t86 * t352) * t237) * t143;
t18 = (t359 + (-t149 + t357) * t66 + (0.2e1 * t48 - t365) * t219) * t96 * t145;
t17 = (t360 + (-t148 + t355) * t65 + (0.2e1 * t47 - t366) * t217) * t95 * t144;
t16 = (t361 + (-t147 + t353) * t64 + (0.2e1 * t46 - t367) * t215) * t94 * t143;
t15 = (-t66 * t285 + (t427 * t66 - t45) * t219 * t225 + (-t201 * t334 + (qJ(3,1) * t313 - t240 + t318) * t225 - t166 * t202 * t310) * t96) * t289 + (-t96 * t285 + (t23 + t286) * t321 + t23 * t166) * t290 + (t166 * t66 - t374) * t291;
t14 = (-t65 * t284 + (t427 * t65 - t44) * t217 * t223 + (-t199 * t335 + (qJ(3,2) * t314 - t238 + t318) * t223 - t164 * t200 * t310) * t95) * t292 + (-t95 * t284 + (t22 + t287) * t322 + t22 * t164) * t293 + (t164 * t65 - t378) * t294;
t13 = (-t64 * t283 + (t427 * t64 - t43) * t215 * t221 + (-t197 * t336 + (qJ(3,3) * t315 - t236 + t318) * t221 - t162 * t198 * t310) * t94) * t295 + (-t94 * t283 + (t24 + t288) * t323 + t24 * t162) * t296 + (t162 * t64 - t382) * t297;
t12 = (-(t205 + t406) * t201 * t362 + (-0.3e1 * (-t240 / 0.3e1 + t205) * t358 + (t240 - t205) * t317) * t423 + ((-t36 * t240 + (-0.4e1 * t306 - t45 + t420) * t235) * t219 + t362 * t406 - t96 * t337) * t225) * t289 + ((t21 * t235 + t286 * t334) * t225 + t21 * pkin(1)) * t290 + (pkin(1) * t66 - t374) * t235 * t291 + (((t45 - 0.2e1 * t420) * t202 + (-t152 - t240) * t383 - 0.2e1 * t306 - t45 + t66 * t426) * t289 + (t21 * t219 + (t409 + 0.1e1) * t96 * t427) * t290 + (t66 * t324 + t402 * t359) * t291) * qJ(3,1);
t11 = (-(t205 + t407) * t199 * t363 + (-0.3e1 * (-t238 / 0.3e1 + t205) * t356 + (t238 - t205) * t317) * t424 + ((-t35 * t238 + (-0.4e1 * t305 - t44 + t421) * t235) * t217 + t363 * t407 - t95 * t337) * t223) * t292 + ((t20 * t235 + t287 * t335) * t223 + t20 * pkin(1)) * t293 + (pkin(1) * t65 - t378) * t235 * t294 + (((t44 - 0.2e1 * t421) * t200 + (-t152 - t238) * t384 - 0.2e1 * t305 - t44 + t65 * t426) * t292 + (t20 * t217 + (t410 + 0.1e1) * t95 * t427) * t293 + (t65 * t325 + t403 * t360) * t294) * qJ(3,2);
t10 = (-(t205 + t408) * t197 * t364 + (-0.3e1 * (-t236 / 0.3e1 + t205) * t354 + (t236 - t205) * t317) * t425 + ((-t34 * t236 + (-0.4e1 * t304 - t43 + t422) * t235) * t215 + t364 * t408 - t94 * t337) * t221) * t295 + ((t19 * t235 + t288 * t336) * t221 + t19 * pkin(1)) * t296 + (pkin(1) * t64 - t382) * t235 * t297 + (((t43 - 0.2e1 * t422) * t198 + (-t152 - t236) * t385 - 0.2e1 * t304 - t43 + t64 * t426) * t295 + (t19 * t215 + (t411 + 0.1e1) * t94 * t427) * t296 + (t64 * t326 + t404 * t361) * t297) * qJ(3,3);
t9 = t15 * t396 - t18 * t298 + (-t12 + t63 * t208 + (-t402 * t208 + (-pkin(1) - t328) * t219) * t93 + g(3) * t225 - t264 * t219) * m(3);
t8 = t14 * t396 - t17 * t299 + (-t11 + t62 * t207 + (-t403 * t207 + (-pkin(1) - t329) * t217) * t92 + g(3) * t223 - t265 * t217) * m(3);
t7 = t13 * t396 - t16 * t300 + (-t10 + t61 * t206 + (-t404 * t206 + (-pkin(1) - t330) * t215) * t91 + g(3) * t221 - t266 * t215) * m(3);
t6 = -t102 * t18 - t105 * t15 + (t264 * t155 - t392) * t225 + (g(3) * t155 + t264 * t159) * t219 + (t228 * t12 - 0.2e1 * t66 * t386) * m(3) + (t136 * t409 + (t108 * t219 + t393) * t225 + t159 * t389 + t136) * t93;
t5 = -t101 * t17 - t104 * t14 + (t265 * t154 - t392) * t223 + (g(3) * t154 + t265 * t159) * t217 + (t228 * t11 - 0.2e1 * t65 * t387) * m(3) + (t135 * t410 + (t107 * t217 + t394) * t223 + t159 * t390 + t135) * t92;
t4 = -t100 * t16 - t103 * t13 + (t266 * t153 - t392) * t221 + (g(3) * t153 + t266 * t159) * t215 + (t228 * t10 - 0.2e1 * t64 * t388) * m(3) + (t134 * t411 + (t106 * t215 + t395) * t221 + t159 * t391 + t134) * t91;
t3 = -t81 * t18 - t102 * t15 - t12 * t298 - 0.4e1 * (-t60 - t307 / 0.2e1) * t423 + (-0.2e1 * (t108 * t66 - t48 * t396) * t383 - t66 * (t133 * t66 + t48 * t316 + t393 * t412)) * t225 + ((t278 + t301) * t63 + (m(3) * t48 - t159 * t66) * t96 * t317) * t219 + (-t60 - t307) * t412 + (t257 * t139 + t142 * t146) * t226 + t220 * (t139 * t146 - t257 * t142);
t2 = -t80 * t17 - t101 * t14 - t11 * t299 - 0.4e1 * (-t59 - t308 / 0.2e1) * t424 + (-0.2e1 * (t107 * t65 - t47 * t396) * t384 - t65 * (t133 * t65 + t47 * t316 + t394 * t413)) * t223 + ((t278 + t302) * t62 + (m(3) * t47 - t159 * t65) * t95 * t317) * t217 + (-t59 - t308) * t413 + (t258 * t138 + t141 * t146) * t224 + t218 * (t138 * t146 - t258 * t141);
t1 = -t79 * t16 - t100 * t13 - t10 * t300 - 0.4e1 * (-t58 - t309 / 0.2e1) * t425 + (-0.2e1 * (t106 * t64 - t46 * t396) * t385 - t64 * (t133 * t64 + t46 * t316 + t395 * t414)) * t221 + ((t278 + t303) * t61 + (m(3) * t46 - t159 * t64) * t94 * t317) * t215 + (-t58 - t309) * t414 + (t259 * t137 + t140 * t146) * t222 + t216 * (t137 * t146 - t259 * t140);
t31 = [(-g(1) + t214) * m(4) + ((t149 * t54 + t219 * t42) * t241 + (t148 * t52 + t217 * t40) * t239 + (t147 * t50 + t215 * t38) * t237) * t212 + (-t30 * t342 + (-t30 * t214 - t3) * t123 + ((-t42 * t373 + t54 * t77) * t214 + (t42 * t372 + t54 * t78) * t213 - t6 * t373 + t77 * t9) * t241) * t145 + (-t28 * t344 + (-t28 * t214 - t2) * t121 + ((-t40 * t377 + t52 * t75) * t214 + (t40 * t376 + t52 * t76) * t213 - t5 * t377 + t75 * t8) * t239) * t144 + (-t26 * t346 + (-t26 * t214 - t1) * t119 + ((-t38 * t381 + t50 * t73) * t214 + (t38 * t380 + t50 * t74) * t213 - t4 * t381 + t73 * t7) * t237) * t143; (-g(2) + t213) * m(4) + ((t149 * t53 + t219 * t41) * t241 + (t148 * t51 + t217 * t39) * t239 + (t147 * t49 + t215 * t37) * t237) * t212 + (-t29 * t341 + (-t29 * t213 - t3) * t122 + ((-t41 * t373 + t53 * t77) * t214 + (t41 * t372 + t53 * t78) * t213 + t6 * t372 + t78 * t9) * t241) * t145 + (-t27 * t343 + (-t27 * t213 - t2) * t120 + ((-t39 * t377 + t51 * t75) * t214 + (t39 * t376 + t51 * t76) * t213 + t5 * t376 + t76 * t8) * t239) * t144 + (-t25 * t345 + (-t25 * t213 - t1) * t118 + ((-t37 * t381 + t49 * t73) * t214 + (t37 * t380 + t49 * t74) * t213 + t4 * t380 + t74 * t7) * t237) * t143; -m(4) * g(3) + (t149 * t9 + t219 * t6) * t241 + (t148 * t8 + t217 * t5) * t239 + (t147 * t7 + t215 * t4) * t237 + (m(4) + (t111 * t149 + t219 * t99) * t241 + (t110 * t148 + t217 * t98) * t239 + (t109 * t147 + t215 * t97) * t237) * t212 + ((-t341 - t342) * (t149 * t397 + t102) * t219 + (t111 * t77 - t89 * t371) * t214 + (t111 * t78 + t90 * t371) * t213) * t338 + ((-t343 - t344) * (t148 * t397 + t101) * t217 + (t110 * t75 - t87 * t375) * t214 + (t110 * t76 + t88 * t375) * t213) * t339 + ((-t345 - t346) * (t147 * t397 + t100) * t215 + (t109 * t73 - t85 * t379) * t214 + (t109 * t74 + t86 * t379) * t213) * t340;];
tauX  = t31;
