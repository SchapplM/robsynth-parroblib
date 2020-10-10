% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G3A0
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:03
% EndTime: 2020-08-06 18:04:13
% DurationCPUTime: 9.10s
% Computational Cost: add. (39630->500), mult. (88452->926), div. (5112->7), fcn. (84471->22), ass. (0->335)
t198 = sin(qJ(2,1));
t204 = cos(qJ(2,1));
t211 = pkin(7) + pkin(6);
t133 = pkin(2) * t198 - t204 * t211;
t184 = sin(pkin(4));
t186 = cos(pkin(4));
t197 = sin(qJ(3,1));
t284 = t186 * t197;
t105 = pkin(3) * t284 + t184 * t133;
t203 = cos(qJ(3,1));
t296 = t184 * t198;
t182 = t203 ^ 2;
t362 = pkin(3) * t182;
t87 = 0.1e1 / (pkin(2) * t284 + t105 * t203 + t296 * t362);
t196 = sin(qJ(2,2));
t202 = cos(qJ(2,2));
t132 = pkin(2) * t196 - t202 * t211;
t195 = sin(qJ(3,2));
t286 = t186 * t195;
t104 = pkin(3) * t286 + t184 * t132;
t201 = cos(qJ(3,2));
t298 = t184 * t196;
t181 = t201 ^ 2;
t363 = pkin(3) * t181;
t86 = 0.1e1 / (pkin(2) * t286 + t104 * t201 + t298 * t363);
t194 = sin(qJ(2,3));
t200 = cos(qJ(2,3));
t131 = pkin(2) * t194 - t200 * t211;
t193 = sin(qJ(3,3));
t288 = t186 * t193;
t103 = pkin(3) * t288 + t184 * t131;
t199 = cos(qJ(3,3));
t300 = t184 * t194;
t180 = t199 ^ 2;
t364 = pkin(3) * t180;
t85 = 0.1e1 / (pkin(2) * t288 + t103 * t199 + t300 * t364);
t189 = legFrame(1,2);
t170 = sin(t189);
t173 = cos(t189);
t126 = t170 * g(1) + t173 * g(2);
t129 = t173 * g(1) - t170 * g(2);
t185 = cos(pkin(8));
t167 = g(3) * t185;
t183 = sin(pkin(8));
t243 = -t129 * t183 - t167;
t396 = t126 * t184 + t243 * t186;
t188 = legFrame(2,2);
t169 = sin(t188);
t172 = cos(t188);
t125 = t169 * g(1) + t172 * g(2);
t128 = t172 * g(1) - t169 * g(2);
t245 = -t128 * t183 - t167;
t395 = t125 * t184 + t245 * t186;
t187 = legFrame(3,2);
t168 = sin(t187);
t171 = cos(t187);
t124 = t168 * g(1) + t171 * g(2);
t127 = t171 * g(1) - t168 * g(2);
t247 = -t127 * t183 - t167;
t394 = t124 * t184 + t247 * t186;
t387 = -rSges(3,1) * t203 + rSges(3,2) * t197;
t386 = -rSges(3,1) * t201 + rSges(3,2) * t195;
t385 = -rSges(3,1) * t199 + rSges(3,2) * t193;
t379 = m(3) * rSges(3,1);
t269 = rSges(3,2) * t379;
t162 = -Icges(3,4) + t269;
t212 = pkin(2) * m(3);
t265 = t212 / 0.2e1;
t252 = rSges(3,1) * t265;
t384 = t162 * t180 + t193 * t252;
t383 = t162 * t181 + t195 * t252;
t382 = t162 * t182 + t197 * t252;
t381 = 0.2e1 * pkin(2);
t208 = xDP(3);
t209 = xDP(2);
t210 = xDP(1);
t238 = t168 * t209 - t171 * t210;
t287 = t186 * t194;
t118 = t183 * t287 - t185 * t200;
t295 = t184 * t199;
t94 = t193 * t118 + t183 * t295;
t121 = t183 * t200 + t185 * t287;
t97 = t193 * t121 + t185 * t295;
t64 = (t208 * t94 + t238 * t97) * t85;
t61 = t64 ^ 2;
t237 = t169 * t209 - t172 * t210;
t285 = t186 * t196;
t119 = t183 * t285 - t185 * t202;
t293 = t184 * t201;
t95 = t195 * t119 + t183 * t293;
t122 = t183 * t202 + t185 * t285;
t98 = t195 * t122 + t185 * t293;
t65 = (t208 * t95 + t237 * t98) * t86;
t62 = t65 ^ 2;
t236 = t170 * t209 - t173 * t210;
t283 = t186 * t198;
t120 = t183 * t283 - t185 * t204;
t291 = t184 * t203;
t96 = t197 * t120 + t183 * t291;
t123 = t183 * t204 + t185 * t283;
t99 = t197 * t123 + t185 * t291;
t66 = (t208 * t96 + t236 * t99) * t87;
t63 = t66 ^ 2;
t380 = -0.2e1 * t162;
t216 = 0.1e1 / pkin(3);
t134 = pkin(2) * t200 + t194 * t211;
t282 = t186 * t200;
t289 = t185 * t186;
t359 = t199 * pkin(3);
t79 = (t183 * t194 - t185 * t282) * t359 - t134 * t289 + t131 * t183;
t302 = t183 * t186;
t82 = (t183 * t282 + t185 * t194) * t359 + t134 * t302 + t131 * t185;
t58 = (t208 * t82 - t238 * t79) * t85 * t216;
t378 = pkin(3) * t58;
t135 = pkin(2) * t202 + t196 * t211;
t281 = t186 * t202;
t357 = t201 * pkin(3);
t80 = (t183 * t196 - t185 * t281) * t357 - t135 * t289 + t132 * t183;
t83 = (t183 * t281 + t185 * t196) * t357 + t135 * t302 + t132 * t185;
t59 = (t208 * t83 - t237 * t80) * t86 * t216;
t377 = pkin(3) * t59;
t136 = pkin(2) * t204 + t198 * t211;
t280 = t186 * t204;
t356 = t203 * pkin(3);
t81 = (t183 * t198 - t185 * t280) * t356 - t136 * t289 + t133 * t183;
t84 = (t183 * t280 + t185 * t198) * t356 + t136 * t302 + t133 * t185;
t60 = (t208 * t84 - t236 * t81) * t87 * t216;
t376 = pkin(3) * t60;
t213 = rSges(3,2) ^ 2;
t214 = rSges(3,1) ^ 2;
t146 = (-t213 + t214) * m(3) + Icges(3,2) - Icges(3,1);
t375 = t146 / 0.2e1;
t205 = pkin(6) + rSges(3,3);
t355 = t205 * m(3);
t150 = rSges(3,2) * t355 - Icges(3,6);
t374 = -t150 / 0.4e1;
t151 = rSges(3,1) * t355 - Icges(3,5);
t373 = t151 / 0.4e1;
t372 = -t162 / 0.2e1;
t335 = t199 * rSges(3,2);
t143 = t193 * rSges(3,1) + t335;
t100 = -t143 * t300 - t186 * t385;
t371 = m(3) * t100;
t333 = t201 * rSges(3,2);
t144 = t195 * rSges(3,1) + t333;
t101 = -t144 * t298 - t186 * t386;
t370 = m(3) * t101;
t331 = t203 * rSges(3,2);
t145 = t197 * rSges(3,1) + t331;
t102 = -t145 * t296 - t186 * t387;
t369 = m(3) * t102;
t368 = m(3) * t216;
t367 = pkin(2) * t193;
t366 = pkin(2) * t195;
t365 = pkin(2) * t197;
t361 = t183 * g(3);
t217 = pkin(2) ^ 2;
t164 = t211 ^ 2 + t217;
t215 = pkin(3) ^ 2;
t338 = t193 * t58;
t268 = pkin(3) * t338;
t275 = pkin(3) * t381;
t360 = (-t211 * t268 + (t180 * t215 + t199 * t275 + t164) * t64) * t64;
t337 = t195 * t59;
t267 = pkin(3) * t337;
t358 = (-t211 * t267 + (t181 * t215 + t201 * t275 + t164) * t65) * t65;
t336 = t197 * t60;
t266 = pkin(3) * t336;
t354 = (-t211 * t266 + (t182 * t215 + t203 * t275 + t164) * t66) * t66;
t347 = t168 * t97;
t346 = t169 * t98;
t345 = t170 * t99;
t344 = t171 * t97;
t343 = t172 * t98;
t342 = t173 * t99;
t158 = m(2) * rSges(2,1) + t212;
t115 = -m(3) * t385 + t158;
t148 = m(2) * rSges(2,2) - t355;
t91 = t115 * t200 - t194 * t148;
t341 = t184 * t91;
t116 = -m(3) * t386 + t158;
t92 = t116 * t202 - t196 * t148;
t340 = t184 * t92;
t117 = -m(3) * t387 + t158;
t93 = t117 * t204 - t198 * t148;
t339 = t184 * t93;
t334 = t200 * t64;
t332 = t202 * t65;
t330 = t204 * t66;
t329 = t216 * t79;
t328 = t216 * t80;
t327 = t216 * t81;
t326 = t216 * t82;
t325 = t216 * t83;
t324 = t216 * t84;
t323 = t64 * t211;
t322 = t65 * t211;
t321 = t66 * t211;
t106 = -t150 * t199 - t151 * t193;
t320 = t106 * t216;
t107 = -t150 * t201 - t151 * t195;
t319 = t107 * t216;
t108 = -t150 * t203 - t151 * t197;
t318 = t108 * t216;
t316 = t124 * t186;
t314 = t125 * t186;
t312 = t126 * t186;
t147 = (t213 + t214) * m(3) + Icges(3,3);
t308 = t147 * t216;
t304 = rSges(3,2) * t184 * t167;
t303 = t183 * t184;
t301 = t184 * t193;
t299 = t184 * t195;
t297 = t184 * t197;
t294 = t184 * t200;
t292 = t184 * t202;
t290 = t184 * t204;
t279 = t186 * t216;
t278 = t194 * t199;
t277 = t196 * t201;
t276 = t198 * t203;
t271 = -t269 / 0.2e1 + Icges(3,4) / 0.2e1;
t270 = -0.2e1 * rSges(3,2) * pkin(2);
t264 = t100 * t368;
t263 = t101 * t368;
t262 = t102 * t368;
t261 = t168 * t329;
t260 = t169 * t328;
t259 = t170 * t327;
t258 = t171 * t329;
t257 = t172 * t328;
t256 = t173 * t327;
t255 = t193 * t323;
t254 = t195 * t322;
t253 = t197 * t321;
t166 = rSges(3,2) * t265;
t251 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t246 = -t127 * t185 + t361;
t244 = -t128 * t185 + t361;
t242 = -t129 * t185 + t361;
t22 = t255 - t378;
t10 = (((t186 * t58 + t64 * t294) * t364 + ((-t268 + t323) * t194 + pkin(2) * t334) * t295 + t22 * t186) * t64 + (t58 * t294 + (t180 * t186 - t278 * t301 - t186) * t64) * t378) * t85;
t112 = pkin(3) * t278 + t131;
t13 = t85 * t279 * t360 + (-t186 * t255 + (-t112 * t301 + (pkin(2) * t199 + t364) * t186) * t58) / (t112 * t295 + (pkin(2) + t359) * t288) * t58;
t16 = (-t199 * t360 - (pkin(2) * t58 - t22 * t199) * t378) * t85;
t149 = t205 ^ 2 + t213 + t217;
t175 = t379 * t381;
t76 = t146 * t180 + (t193 * t380 + t175) * t199 + (t193 * t270 + t149) * m(3) + t251;
t4 = -t16 * t341 - t76 * t10 - t106 * t13 - 0.4e1 * ((t374 * t193 + t373 * t199) * t58 + ((t193 * t375 + t166) * t199 + t372 + t384) * t64) * t58 + (-t115 * t394 - t246 * t148) * t200 - (t246 * t115 - t148 * t394) * t194;
t130 = (t214 / 0.2e1 - t213 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t220 = t394 * t194 - t246 * t200;
t7 = -t16 * t371 - t106 * t10 - t147 * t13 + 0.2e1 * ((t130 * t193 + t166) * t199 + t271 + t384) * t61 + m(3) * (((t247 * t184 - t316) * rSges(3,1) + t220 * rSges(3,2)) * t199 + (t304 + (t127 * t303 + t316) * rSges(3,2) + t220 * rSges(3,1)) * t193);
t235 = t7 * t329 - t97 * t4;
t23 = t254 - t377;
t11 = (((t186 * t59 + t65 * t292) * t363 + ((-t267 + t322) * t196 + pkin(2) * t332) * t293 + t23 * t186) * t65 + (t59 * t292 + (t181 * t186 - t277 * t299 - t186) * t65) * t377) * t86;
t113 = pkin(3) * t277 + t132;
t14 = t86 * t279 * t358 + (-t186 * t254 + (-t113 * t299 + (pkin(2) * t201 + t363) * t186) * t59) / (t113 * t293 + (pkin(2) + t357) * t286) * t59;
t17 = (-t201 * t358 - (pkin(2) * t59 - t23 * t201) * t377) * t86;
t77 = t146 * t181 + (t195 * t380 + t175) * t201 + (t195 * t270 + t149) * m(3) + t251;
t5 = -t17 * t340 - t77 * t11 - t107 * t14 - 0.4e1 * ((t374 * t195 + t373 * t201) * t59 + ((t195 * t375 + t166) * t201 + t372 + t383) * t65) * t59 + (-t116 * t395 - t244 * t148) * t202 - (t244 * t116 - t148 * t395) * t196;
t219 = t395 * t196 - t244 * t202;
t8 = -t17 * t370 - t107 * t11 - t147 * t14 + 0.2e1 * ((t130 * t195 + t166) * t201 + t271 + t383) * t62 + m(3) * (((t245 * t184 - t314) * rSges(3,1) + t219 * rSges(3,2)) * t201 + (t304 + (t128 * t303 + t314) * rSges(3,2) + t219 * rSges(3,1)) * t195);
t234 = t8 * t328 - t98 * t5;
t24 = t253 - t376;
t12 = (((t186 * t60 + t66 * t290) * t362 + ((-t266 + t321) * t198 + pkin(2) * t330) * t291 + t24 * t186) * t66 + (t60 * t290 + (t182 * t186 - t276 * t297 - t186) * t66) * t376) * t87;
t114 = pkin(3) * t276 + t133;
t15 = t87 * t279 * t354 + (-t186 * t253 + (-t114 * t297 + (pkin(2) * t203 + t362) * t186) * t60) / (t114 * t291 + (pkin(2) + t356) * t284) * t60;
t18 = (-t203 * t354 - (pkin(2) * t60 - t24 * t203) * t376) * t87;
t78 = t146 * t182 + (t197 * t380 + t175) * t203 + (t197 * t270 + t149) * m(3) + t251;
t6 = -t18 * t339 - t78 * t12 - t108 * t15 - 0.4e1 * ((t374 * t197 + t373 * t203) * t60 + ((t197 * t375 + t166) * t203 + t372 + t382) * t66) * t60 + (-t117 * t396 - t242 * t148) * t204 - (t242 * t117 - t148 * t396) * t198;
t218 = t396 * t198 - t242 * t204;
t9 = -t18 * t369 - t108 * t12 - t147 * t15 + 0.2e1 * ((t130 * t197 + t166) * t203 + t271 + t382) * t63 + m(3) * (((t243 * t184 - t312) * rSges(3,1) + t218 * rSges(3,2)) * t203 + (t304 + (t129 * t303 + t312) * rSges(3,2) + t218 * rSges(3,1)) * t197);
t233 = t9 * t327 - t99 * t6;
t232 = t79 * t320 - t76 * t97;
t231 = t80 * t319 - t77 * t98;
t230 = t81 * t318 - t78 * t99;
t229 = pkin(3) * t301 - t131 * t186;
t228 = pkin(3) * t299 - t132 * t186;
t227 = pkin(3) * t297 - t133 * t186;
t226 = -t106 * t97 + t79 * t308;
t225 = -t107 * t98 + t80 * t308;
t224 = -t108 * t99 + t81 * t308;
t223 = t79 * t264 - t97 * t341;
t222 = t80 * t263 - t98 * t340;
t221 = t81 * t262 - t99 * t339;
t192 = xDDP(1);
t191 = xDDP(2);
t190 = xDDP(3);
t179 = m(1) + m(2) + m(3);
t90 = t185 * t136 + t227 * t183;
t89 = t185 * t135 + t228 * t183;
t88 = t185 * t134 + t229 * t183;
t75 = -t123 * t362 - t136 * t183 * t203 + (pkin(2) * t297 + t227 * t203) * t185;
t74 = -t122 * t363 - t135 * t183 * t201 + (pkin(2) * t299 + t228 * t201) * t185;
t73 = -t121 * t364 - t134 * t183 * t199 + (pkin(2) * t301 + t229 * t199) * t185;
t72 = (t120 * t170 + t173 * t296) * t362 + (t173 * t105 - t90 * t170) * t203 + (-t170 * t303 + t173 * t186) * t365;
t71 = (t119 * t169 + t172 * t298) * t363 + (t172 * t104 - t89 * t169) * t201 + (-t169 * t303 + t172 * t186) * t366;
t70 = (t118 * t168 + t171 * t300) * t364 + (t171 * t103 - t88 * t168) * t199 + (-t168 * t303 + t171 * t186) * t367;
t69 = -(t120 * t173 - t170 * t296) * t362 + (t170 * t105 + t90 * t173) * t203 + (t186 * t170 + t173 * t303) * t365;
t68 = -(t119 * t172 - t169 * t298) * t363 + (t169 * t104 + t89 * t172) * t201 + (t186 * t169 + t172 * t303) * t366;
t67 = -(t118 * t171 - t168 * t300) * t364 + (t168 * t103 + t88 * t171) * t199 + (t186 * t168 + t171 * t303) * t367;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = (t108 * t96 + t84 * t308 + t75 * t369) * t87;
t53 = (t107 * t95 + t83 * t308 + t74 * t370) * t86;
t52 = (t106 * t94 + t82 * t308 + t73 * t371) * t85;
t51 = (t179 * t75 + t84 * t262 + t96 * t339) * t87;
t50 = (t179 * t74 + t83 * t263 + t95 * t340) * t86;
t49 = (t179 * t73 + t82 * t264 + t94 * t341) * t85;
t48 = (t84 * t318 + t75 * t339 + t78 * t96) * t87;
t47 = (t83 * t319 + t74 * t340 + t77 * t95) * t86;
t46 = (t82 * t320 + t73 * t341 + t76 * t94) * t85;
t45 = (-t224 * t170 + t72 * t369) * t87;
t44 = (-t225 * t169 + t71 * t370) * t86;
t43 = (-t226 * t168 + t70 * t371) * t85;
t42 = (t224 * t173 + t69 * t369) * t87;
t41 = (t225 * t172 + t68 * t370) * t86;
t40 = (t226 * t171 + t67 * t371) * t85;
t39 = (-t221 * t170 + t179 * t72) * t87;
t38 = (-t222 * t169 + t179 * t71) * t86;
t37 = (-t223 * t168 + t179 * t70) * t85;
t36 = (t221 * t173 + t179 * t69) * t87;
t35 = (t222 * t172 + t179 * t68) * t86;
t34 = (t223 * t171 + t179 * t67) * t85;
t33 = (-t230 * t170 + t72 * t339) * t87;
t32 = (-t231 * t169 + t71 * t340) * t86;
t31 = (-t232 * t168 + t70 * t341) * t85;
t30 = (t230 * t173 + t69 * t339) * t87;
t29 = (t231 * t172 + t68 * t340) * t86;
t28 = (t232 * t171 + t67 * t341) * t85;
t3 = (-t93 * t12 + (-t148 * t204 - t158 * t198) * t63) * t184 + (-t18 - t126) * t179 + (-t102 * t15 + (-0.2e1 * (rSges(3,1) * t336 + t60 * t331) * t330 + t387 * t198 * (t63 + t57)) * t184 - t57 * t186 * t145) * m(3);
t2 = (-t92 * t11 + (-t148 * t202 - t158 * t196) * t62) * t184 + (-t17 - t125) * t179 + (-t101 * t14 + (-0.2e1 * (rSges(3,1) * t337 + t59 * t333) * t332 + t386 * t196 * (t62 + t56)) * t184 - t56 * t186 * t144) * m(3);
t1 = (-t91 * t10 + (-t148 * t200 - t158 * t194) * t61) * t184 + (-t16 - t124) * t179 + (-t100 * t13 + (-0.2e1 * (rSges(3,1) * t338 + t58 * t335) * t334 + t385 * t194 * (t61 + t55)) * t184 - t55 * t186 * t143) * m(3);
t19 = [(-g(1) + t192) * m(4) + ((-t42 * t259 + t30 * t345 + t36 * t72) * t191 + (t30 * t96 + t42 * t324 + t36 * t75) * t190 + (t36 * t192 + t3) * t69 + ((-t30 * t99 + t42 * t327) * t192 + t233) * t173) * t87 + ((-t41 * t260 + t29 * t346 + t35 * t71) * t191 + (t29 * t95 + t41 * t325 + t35 * t74) * t190 + (t35 * t192 + t2) * t68 + ((-t29 * t98 + t41 * t328) * t192 + t234) * t172) * t86 + ((-t40 * t261 + t28 * t347 + t34 * t70) * t191 + (t28 * t94 + t40 * t326 + t34 * t73) * t190 + (t34 * t192 + t1) * t67 + ((-t28 * t97 + t40 * t329) * t192 + t235) * t171) * t85; (-g(2) + t191) * m(4) + ((t45 * t256 - t33 * t342 + t39 * t69) * t192 + (t45 * t324 + t33 * t96 + t39 * t75) * t190 + (t39 * t191 + t3) * t72 + ((-t45 * t327 + t33 * t99) * t191 - t233) * t170) * t87 + ((t44 * t257 - t32 * t343 + t38 * t68) * t192 + (t32 * t95 + t44 * t325 + t38 * t74) * t190 + (t38 * t191 + t2) * t71 + ((t32 * t98 - t44 * t328) * t191 - t234) * t169) * t86 + ((t43 * t258 - t31 * t344 + t37 * t67) * t192 + (t31 * t94 + t43 * t326 + t37 * t73) * t190 + (t37 * t191 + t1) * t70 + ((t31 * t97 - t43 * t329) * t191 - t235) * t168) * t85; (-g(3) + t190) * m(4) + ((t54 * t256 - t48 * t342 + t51 * t69) * t192 + (-t54 * t259 + t48 * t345 + t51 * t72) * t191 + (t54 * t324 + t48 * t96 + t51 * t75) * t190 + t75 * t3 + t96 * t6 + t9 * t324) * t87 + ((t53 * t257 - t47 * t343 + t50 * t68) * t192 + (-t53 * t260 + t47 * t346 + t50 * t71) * t191 + (t53 * t325 + t47 * t95 + t50 * t74) * t190 + t74 * t2 + t95 * t5 + t8 * t325) * t86 + ((t52 * t258 - t46 * t344 + t49 * t67) * t192 + (-t52 * t261 + t46 * t347 + t49 * t70) * t191 + (t52 * t326 + t46 * t94 + t49 * t73) * t190 + t73 * t1 + t94 * t4 + t7 * t326) * t85;];
tauX  = t19;
