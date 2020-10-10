% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G2A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:48:58
% EndTime: 2020-08-06 17:49:08
% DurationCPUTime: 9.09s
% Computational Cost: add. (39630->503), mult. (88452->927), div. (5112->7), fcn. (84471->22), ass. (0->342)
t199 = sin(qJ(2,1));
t205 = cos(qJ(2,1));
t212 = pkin(7) + pkin(6);
t134 = pkin(2) * t199 - t205 * t212;
t185 = sin(pkin(4));
t187 = cos(pkin(4));
t198 = sin(qJ(3,1));
t278 = t198 * t187;
t105 = pkin(3) * t278 + t134 * t185;
t204 = cos(qJ(3,1));
t300 = t185 * t199;
t183 = t204 ^ 2;
t368 = pkin(3) * t183;
t87 = 0.1e1 / (pkin(2) * t278 + t105 * t204 + t300 * t368);
t197 = sin(qJ(2,2));
t203 = cos(qJ(2,2));
t133 = pkin(2) * t197 - t203 * t212;
t196 = sin(qJ(3,2));
t280 = t196 * t187;
t104 = pkin(3) * t280 + t133 * t185;
t202 = cos(qJ(3,2));
t302 = t185 * t197;
t182 = t202 ^ 2;
t369 = pkin(3) * t182;
t86 = 0.1e1 / (pkin(2) * t280 + t104 * t202 + t302 * t369);
t195 = sin(qJ(2,3));
t201 = cos(qJ(2,3));
t132 = pkin(2) * t195 - t201 * t212;
t194 = sin(qJ(3,3));
t282 = t194 * t187;
t103 = pkin(3) * t282 + t132 * t185;
t200 = cos(qJ(3,3));
t304 = t185 * t195;
t181 = t200 ^ 2;
t370 = pkin(3) * t181;
t85 = 0.1e1 / (pkin(2) * t282 + t103 * t200 + t304 * t370);
t190 = legFrame(1,2);
t171 = sin(t190);
t174 = cos(t190);
t126 = g(1) * t171 + g(2) * t174;
t129 = g(1) * t174 - g(2) * t171;
t184 = sin(pkin(8));
t168 = g(3) * t184;
t186 = cos(pkin(8));
t243 = t129 * t186 - t168;
t403 = t126 * t185 + t243 * t187;
t189 = legFrame(2,2);
t170 = sin(t189);
t173 = cos(t189);
t125 = g(1) * t170 + g(2) * t173;
t128 = g(1) * t173 - g(2) * t170;
t245 = t128 * t186 - t168;
t402 = t125 * t185 + t245 * t187;
t188 = legFrame(3,2);
t169 = sin(t188);
t172 = cos(t188);
t124 = g(1) * t169 + g(2) * t172;
t127 = g(1) * t172 - g(2) * t169;
t247 = t127 * t186 - t168;
t401 = t124 * t185 + t247 * t187;
t394 = -rSges(3,1) * t204 + rSges(3,2) * t198;
t393 = -rSges(3,1) * t202 + rSges(3,2) * t196;
t392 = -rSges(3,1) * t200 + rSges(3,2) * t194;
t386 = m(3) * rSges(3,1);
t270 = rSges(3,2) * t386;
t163 = -Icges(3,4) + t270;
t213 = pkin(2) * m(3);
t266 = t213 / 0.2e1;
t253 = rSges(3,1) * t266;
t391 = t163 * t181 + t194 * t253;
t390 = t163 * t182 + t196 * t253;
t389 = t163 * t183 + t198 * t253;
t388 = 0.2e1 * pkin(2);
t209 = xDP(3);
t210 = xDP(2);
t211 = xDP(1);
t239 = t169 * t210 - t172 * t211;
t289 = t187 * t195;
t118 = t184 * t289 - t186 * t201;
t299 = t185 * t200;
t94 = t118 * t194 + t184 * t299;
t121 = t184 * t201 + t186 * t289;
t292 = t186 * t200;
t97 = -t121 * t194 - t185 * t292;
t64 = (t209 * t97 + t239 * t94) * t85;
t61 = t64 ^ 2;
t238 = t170 * t210 - t173 * t211;
t288 = t187 * t197;
t119 = t184 * t288 - t186 * t203;
t297 = t185 * t202;
t95 = t119 * t196 + t184 * t297;
t122 = t184 * t203 + t186 * t288;
t291 = t186 * t202;
t98 = -t122 * t196 - t185 * t291;
t65 = (t209 * t98 + t238 * t95) * t86;
t62 = t65 ^ 2;
t237 = t171 * t210 - t174 * t211;
t287 = t187 * t199;
t120 = t184 * t287 - t186 * t205;
t295 = t185 * t204;
t96 = t120 * t198 + t184 * t295;
t123 = t184 * t205 + t186 * t287;
t290 = t186 * t204;
t99 = -t123 * t198 - t185 * t290;
t66 = (t209 * t99 + t237 * t96) * t87;
t63 = t66 ^ 2;
t387 = -0.2e1 * t163;
t217 = 0.1e1 / pkin(3);
t135 = pkin(2) * t201 + t212 * t195;
t286 = t187 * t201;
t293 = t186 * t187;
t367 = pkin(3) * t200;
t79 = (t184 * t195 - t186 * t286) * t367 - t135 * t293 + t184 * t132;
t307 = t184 * t187;
t82 = (t184 * t286 + t186 * t195) * t367 + t135 * t307 + t132 * t186;
t58 = (t209 * t79 + t239 * t82) * t85 * t217;
t385 = pkin(3) * t58;
t136 = pkin(2) * t203 + t212 * t197;
t285 = t187 * t203;
t366 = pkin(3) * t202;
t80 = (t184 * t197 - t186 * t285) * t366 - t136 * t293 + t184 * t133;
t83 = (t184 * t285 + t186 * t197) * t366 + t136 * t307 + t133 * t186;
t59 = (t209 * t80 + t238 * t83) * t86 * t217;
t384 = pkin(3) * t59;
t137 = pkin(2) * t205 + t212 * t199;
t284 = t187 * t205;
t365 = pkin(3) * t204;
t81 = (t184 * t199 - t186 * t284) * t365 - t137 * t293 + t184 * t134;
t84 = (t184 * t284 + t186 * t199) * t365 + t137 * t307 + t134 * t186;
t60 = (t209 * t81 + t237 * t84) * t87 * t217;
t383 = pkin(3) * t60;
t214 = rSges(3,2) ^ 2;
t215 = rSges(3,1) ^ 2;
t147 = (-t214 + t215) * m(3) + Icges(3,2) - Icges(3,1);
t382 = t147 / 0.2e1;
t206 = pkin(6) + rSges(3,3);
t375 = m(3) * t206;
t151 = rSges(3,2) * t375 - Icges(3,6);
t381 = -t151 / 0.4e1;
t152 = rSges(3,1) * t375 - Icges(3,5);
t380 = t152 / 0.4e1;
t379 = -t163 / 0.2e1;
t354 = rSges(3,2) * t200;
t144 = rSges(3,1) * t194 + t354;
t100 = -t144 * t304 - t187 * t392;
t378 = m(3) * t100;
t353 = rSges(3,2) * t202;
t145 = rSges(3,1) * t196 + t353;
t101 = -t145 * t302 - t187 * t393;
t377 = m(3) * t101;
t352 = rSges(3,2) * t204;
t146 = rSges(3,1) * t198 + t352;
t102 = -t146 * t300 - t187 * t394;
t376 = m(3) * t102;
t374 = m(3) * t217;
t373 = pkin(2) * t194;
t372 = pkin(2) * t196;
t371 = pkin(2) * t198;
t364 = g(3) * t186;
t218 = pkin(2) ^ 2;
t165 = t212 ^ 2 + t218;
t216 = pkin(3) ^ 2;
t342 = t194 * t58;
t269 = pkin(3) * t342;
t276 = pkin(3) * t388;
t363 = (-t212 * t269 + (t181 * t216 + t200 * t276 + t165) * t64) * t64;
t341 = t196 * t59;
t268 = pkin(3) * t341;
t362 = (-t212 * t268 + (t182 * t216 + t202 * t276 + t165) * t65) * t65;
t340 = t198 * t60;
t267 = pkin(3) * t340;
t361 = (-t212 * t267 + (t183 * t216 + t204 * t276 + t165) * t66) * t66;
t351 = t169 * t94;
t350 = t170 * t95;
t349 = t171 * t96;
t348 = t172 * t94;
t347 = t173 * t95;
t346 = t174 * t96;
t159 = m(2) * rSges(2,1) + t213;
t115 = -m(3) * t392 + t159;
t149 = m(2) * rSges(2,2) - t375;
t91 = t115 * t201 - t149 * t195;
t345 = t185 * t91;
t116 = -m(3) * t393 + t159;
t92 = t116 * t203 - t149 * t197;
t344 = t185 * t92;
t117 = -m(3) * t394 + t159;
t93 = t117 * t205 - t149 * t199;
t343 = t185 * t93;
t339 = t212 * t64;
t338 = t212 * t65;
t337 = t212 * t66;
t336 = t217 * t79;
t335 = t217 * t80;
t334 = t217 * t81;
t333 = t217 * t82;
t332 = t217 * t83;
t331 = t217 * t84;
t330 = t64 * t201;
t329 = t65 * t203;
t328 = t66 * t205;
t106 = -t151 * t200 - t152 * t194;
t327 = t106 * t217;
t107 = -t151 * t202 - t152 * t196;
t326 = t107 * t217;
t108 = -t151 * t204 - t152 * t198;
t325 = t108 * t217;
t323 = t124 * t187;
t321 = t125 * t187;
t319 = t126 * t187;
t148 = (t214 + t215) * m(3) + Icges(3,3);
t315 = t148 * t217;
t311 = rSges(3,2) * t168 * t185;
t310 = t184 * t127;
t309 = t184 * t128;
t308 = t184 * t129;
t306 = t185 * t186;
t305 = t185 * t194;
t303 = t185 * t196;
t301 = t185 * t198;
t298 = t185 * t201;
t296 = t185 * t203;
t294 = t185 * t205;
t283 = t187 * t217;
t281 = t195 * t200;
t279 = t197 * t202;
t277 = t199 * t204;
t272 = -t270 / 0.2e1 + Icges(3,4) / 0.2e1;
t271 = -0.2e1 * rSges(3,2) * pkin(2);
t265 = t100 * t374;
t264 = t101 * t374;
t263 = t102 * t374;
t262 = t169 * t333;
t261 = t170 * t332;
t260 = t171 * t331;
t259 = t172 * t333;
t258 = t173 * t332;
t257 = t174 * t331;
t256 = t194 * t339;
t255 = t196 * t338;
t254 = t198 * t337;
t167 = rSges(3,2) * t266;
t252 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t248 = t310 + t364;
t246 = t309 + t364;
t244 = t308 + t364;
t22 = t256 - t385;
t10 = (((t187 * t58 + t64 * t298) * t370 + ((-t269 + t339) * t195 + pkin(2) * t330) * t299 + t22 * t187) * t64 + (t58 * t298 + (t181 * t187 - t281 * t305 - t187) * t64) * t385) * t85;
t112 = pkin(3) * t281 + t132;
t13 = t85 * t283 * t363 + (-t187 * t256 + (-t112 * t305 + t187 * (pkin(2) * t200 + t370)) * t58) / (t112 * t299 + (pkin(2) + t367) * t282) * t58;
t130 = t149 * t364;
t16 = (-t200 * t363 - (pkin(2) * t58 - t200 * t22) * t385) * t85;
t150 = t206 ^ 2 + t214 + t218;
t176 = t386 * t388;
t76 = t147 * t181 + (t194 * t387 + t176) * t200 + (t194 * t271 + t150) * m(3) + t252;
t4 = -t16 * t345 - t76 * t10 - t106 * t13 - 0.4e1 * ((t381 * t194 + t380 * t200) * t58 + ((t194 * t382 + t167) * t200 + t379 + t391) * t64) * t58 + (-t115 * t401 + t149 * t310 + t130) * t201 - (-t248 * t115 - t149 * t401) * t195;
t131 = (t215 / 0.2e1 - t214 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t221 = t401 * t195 + t248 * t201;
t7 = -t16 * t378 - t106 * t10 - t148 * t13 + 0.2e1 * ((t131 * t194 + t167) * t200 + t272 + t391) * t61 + m(3) * (((t247 * t185 - t323) * rSges(3,1) + t221 * rSges(3,2)) * t200 + t194 * (t311 + (-t127 * t306 + t323) * rSges(3,2) + t221 * rSges(3,1)));
t236 = t7 * t333 + t94 * t4;
t23 = t255 - t384;
t11 = (((t187 * t59 + t65 * t296) * t369 + ((-t268 + t338) * t197 + pkin(2) * t329) * t297 + t23 * t187) * t65 + (t59 * t296 + (t182 * t187 - t279 * t303 - t187) * t65) * t384) * t86;
t113 = pkin(3) * t279 + t133;
t14 = t86 * t283 * t362 + (-t187 * t255 + (-t113 * t303 + t187 * (pkin(2) * t202 + t369)) * t59) / (t113 * t297 + (pkin(2) + t366) * t280) * t59;
t17 = (-t202 * t362 - (pkin(2) * t59 - t202 * t23) * t384) * t86;
t77 = t147 * t182 + (t196 * t387 + t176) * t202 + (t196 * t271 + t150) * m(3) + t252;
t5 = -t17 * t344 - t77 * t11 - t107 * t14 - 0.4e1 * ((t381 * t196 + t380 * t202) * t59 + ((t196 * t382 + t167) * t202 + t379 + t390) * t65) * t59 + (-t116 * t402 + t149 * t309 + t130) * t203 - (-t246 * t116 - t149 * t402) * t197;
t220 = t402 * t197 + t246 * t203;
t8 = -t17 * t377 - t107 * t11 - t148 * t14 + 0.2e1 * ((t131 * t196 + t167) * t202 + t272 + t390) * t62 + m(3) * (((t245 * t185 - t321) * rSges(3,1) + t220 * rSges(3,2)) * t202 + t196 * (t311 + (-t128 * t306 + t321) * rSges(3,2) + t220 * rSges(3,1)));
t235 = t8 * t332 + t95 * t5;
t24 = t254 - t383;
t12 = (((t187 * t60 + t66 * t294) * t368 + ((-t267 + t337) * t199 + pkin(2) * t328) * t295 + t24 * t187) * t66 + (t60 * t294 + (t183 * t187 - t277 * t301 - t187) * t66) * t383) * t87;
t114 = pkin(3) * t277 + t134;
t15 = t87 * t283 * t361 + (-t187 * t254 + (-t114 * t301 + t187 * (pkin(2) * t204 + t368)) * t60) / (t114 * t295 + (pkin(2) + t365) * t278) * t60;
t18 = (-t204 * t361 - (pkin(2) * t60 - t204 * t24) * t383) * t87;
t78 = t147 * t183 + (t198 * t387 + t176) * t204 + (t198 * t271 + t150) * m(3) + t252;
t6 = -t18 * t343 - t78 * t12 - t108 * t15 - 0.4e1 * ((t381 * t198 + t380 * t204) * t60 + ((t198 * t382 + t167) * t204 + t379 + t389) * t66) * t60 + (-t117 * t403 + t149 * t308 + t130) * t205 - (-t244 * t117 - t149 * t403) * t199;
t219 = t403 * t199 + t244 * t205;
t9 = -t18 * t376 - t108 * t12 - t148 * t15 + 0.2e1 * ((t131 * t198 + t167) * t204 + t272 + t389) * t63 + m(3) * (((t243 * t185 - t319) * rSges(3,1) + t219 * rSges(3,2)) * t204 + t198 * (t311 + (-t129 * t306 + t319) * rSges(3,2) + t219 * rSges(3,1)));
t234 = t9 * t331 + t96 * t6;
t233 = t82 * t327 + t76 * t94;
t232 = t83 * t326 + t77 * t95;
t231 = t84 * t325 + t78 * t96;
t230 = pkin(3) * t305 - t132 * t187;
t229 = pkin(3) * t303 - t133 * t187;
t228 = pkin(3) * t301 - t134 * t187;
t227 = t106 * t94 + t82 * t315;
t226 = t107 * t95 + t83 * t315;
t225 = t108 * t96 + t84 * t315;
t224 = t82 * t265 + t94 * t345;
t223 = t83 * t264 + t95 * t344;
t222 = t84 * t263 + t96 * t343;
t193 = xDDP(1);
t192 = xDDP(2);
t191 = xDDP(3);
t180 = m(1) + m(2) + m(3);
t90 = -t184 * t137 + t228 * t186;
t89 = -t184 * t136 + t229 * t186;
t88 = -t184 * t135 + t230 * t186;
t75 = -t120 * t368 + t137 * t290 + (pkin(2) * t301 + t228 * t204) * t184;
t74 = -t119 * t369 + t136 * t291 + (pkin(2) * t303 + t229 * t202) * t184;
t73 = -t118 * t370 + t135 * t292 + (pkin(2) * t305 + t230 * t200) * t184;
t72 = (t123 * t174 + t171 * t300) * t368 + (t105 * t171 - t174 * t90) * t204 + (t171 * t187 - t174 * t306) * t371;
t71 = -(t123 * t171 - t174 * t300) * t368 + (t105 * t174 + t171 * t90) * t204 + (t171 * t306 + t174 * t187) * t371;
t70 = (t122 * t173 + t170 * t302) * t369 + (t104 * t170 - t173 * t89) * t202 + (t170 * t187 - t173 * t306) * t372;
t69 = -(t122 * t170 - t173 * t302) * t369 + (t104 * t173 + t170 * t89) * t202 + (t170 * t306 + t173 * t187) * t372;
t68 = (t121 * t172 + t169 * t304) * t370 + (t103 * t169 - t172 * t88) * t200 + (t169 * t187 - t172 * t306) * t373;
t67 = -(t121 * t169 - t172 * t304) * t370 + (t103 * t172 + t169 * t88) * t200 + (t169 * t306 + t172 * t187) * t373;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = (t108 * t99 + t81 * t315 + t75 * t376) * t87;
t53 = (t107 * t98 + t80 * t315 + t74 * t377) * t86;
t52 = (t106 * t97 + t79 * t315 + t73 * t378) * t85;
t51 = (t180 * t75 + t81 * t263 + t99 * t343) * t87;
t50 = (t180 * t74 + t80 * t264 + t98 * t344) * t86;
t49 = (t180 * t73 + t79 * t265 + t97 * t345) * t85;
t48 = (t81 * t325 + t75 * t343 + t78 * t99) * t87;
t47 = (t80 * t326 + t74 * t344 + t77 * t98) * t86;
t46 = (t79 * t327 + t73 * t345 + t76 * t97) * t85;
t45 = (-t225 * t174 + t72 * t376) * t87;
t44 = (t225 * t171 + t71 * t376) * t87;
t43 = (-t226 * t173 + t70 * t377) * t86;
t42 = (t226 * t170 + t69 * t377) * t86;
t41 = (-t227 * t172 + t68 * t378) * t85;
t40 = (t227 * t169 + t67 * t378) * t85;
t39 = (-t222 * t174 + t180 * t72) * t87;
t38 = (t222 * t171 + t180 * t71) * t87;
t37 = (-t223 * t173 + t180 * t70) * t86;
t36 = (t223 * t170 + t180 * t69) * t86;
t35 = (-t224 * t172 + t180 * t68) * t85;
t34 = (t224 * t169 + t180 * t67) * t85;
t33 = (-t231 * t174 + t72 * t343) * t87;
t32 = (t231 * t171 + t71 * t343) * t87;
t31 = (-t232 * t173 + t70 * t344) * t86;
t30 = (t232 * t170 + t69 * t344) * t86;
t29 = (-t233 * t172 + t68 * t345) * t85;
t28 = (t233 * t169 + t67 * t345) * t85;
t3 = (-t93 * t12 + (-t149 * t205 - t159 * t199) * t63) * t185 + (-t18 - t126) * t180 + (-t102 * t15 + (-0.2e1 * (rSges(3,1) * t340 + t60 * t352) * t328 + t394 * t199 * (t63 + t57)) * t185 - t57 * t187 * t146) * m(3);
t2 = (-t92 * t11 + (-t149 * t203 - t159 * t197) * t62) * t185 + (-t17 - t125) * t180 + (-t101 * t14 + (-0.2e1 * (rSges(3,1) * t341 + t59 * t353) * t329 + t393 * t197 * (t62 + t56)) * t185 - t56 * t187 * t145) * m(3);
t1 = (-t91 * t10 + (-t149 * t201 - t159 * t195) * t61) * t185 + (-t16 - t124) * t180 + (-t100 * t13 + (-0.2e1 * (rSges(3,1) * t342 + t58 * t354) * t330 + t392 * t195 * (t61 + t55)) * t185 - t55 * t187 * t144) * m(3);
t19 = [(-g(1) + t193) * m(4) + ((t45 * t260 + t33 * t349 + t39 * t71) * t192 + (t33 * t99 + t45 * t334 + t39 * t75) * t191 + (t39 * t193 + t3) * t72 + ((-t33 * t96 - t45 * t331) * t193 - t234) * t174) * t87 + ((t43 * t261 + t31 * t350 + t37 * t69) * t192 + (t31 * t98 + t43 * t335 + t37 * t74) * t191 + (t37 * t193 + t2) * t70 + ((-t31 * t95 - t43 * t332) * t193 - t235) * t173) * t86 + ((t41 * t262 + t29 * t351 + t35 * t67) * t192 + (t29 * t97 + t41 * t336 + t35 * t73) * t191 + (t35 * t193 + t1) * t68 + ((-t29 * t94 - t41 * t333) * t193 - t236) * t172) * t85; (-g(2) + t192) * m(4) + ((-t44 * t257 - t32 * t346 + t38 * t72) * t193 + (t32 * t99 + t44 * t334 + t38 * t75) * t191 + (t38 * t192 + t3) * t71 + ((t32 * t96 + t44 * t331) * t192 + t234) * t171) * t87 + ((-t42 * t258 - t30 * t347 + t36 * t70) * t193 + (t30 * t98 + t42 * t335 + t36 * t74) * t191 + (t36 * t192 + t2) * t69 + ((t30 * t95 + t42 * t332) * t192 + t235) * t170) * t86 + ((-t40 * t259 - t28 * t348 + t34 * t68) * t193 + (t28 * t97 + t40 * t336 + t34 * t73) * t191 + (t34 * t192 + t1) * t67 + ((t28 * t94 + t40 * t333) * t192 + t236) * t169) * t85; (-g(3) + t191) * m(4) + ((-t54 * t257 - t48 * t346 + t51 * t72) * t193 + (t54 * t260 + t48 * t349 + t51 * t71) * t192 + (t54 * t334 + t48 * t99 + t51 * t75) * t191 + t75 * t3 + t99 * t6 + t9 * t334) * t87 + ((-t53 * t258 - t47 * t347 + t50 * t70) * t193 + (t53 * t261 + t47 * t350 + t50 * t69) * t192 + (t53 * t335 + t47 * t98 + t50 * t74) * t191 + t74 * t2 + t98 * t5 + t8 * t335) * t86 + ((-t52 * t259 - t46 * t348 + t49 * t68) * t193 + (t52 * t262 + t46 * t351 + t49 * t67) * t192 + (t52 * t336 + t46 * t97 + t49 * t73) * t191 + t73 * t1 + t97 * t4 + t7 * t336) * t85;];
tauX  = t19;
