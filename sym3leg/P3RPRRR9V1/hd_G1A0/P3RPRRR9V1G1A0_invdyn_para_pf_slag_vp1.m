% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR9V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:46:42
% EndTime: 2020-08-06 18:46:49
% DurationCPUTime: 7.02s
% Computational Cost: add. (19797->487), mult. (19257->801), div. (3780->11), fcn. (13953->55), ass. (0->332)
t202 = m(2) + m(3);
t400 = pkin(1) * t202;
t169 = (pkin(7) + qJ(3,1));
t122 = 2 * t169;
t113 = cos(t122);
t399 = (t113 + 0.1e1) * pkin(3);
t168 = (pkin(7) + qJ(3,2));
t121 = 2 * t168;
t112 = cos(t121);
t398 = (t112 + 0.1e1) * pkin(3);
t167 = (pkin(7) + qJ(3,3));
t120 = 2 * t167;
t111 = cos(t120);
t397 = (t111 + 0.1e1) * pkin(3);
t396 = 2 * pkin(1);
t395 = -pkin(3) / 0.2e1;
t209 = 2 * pkin(7);
t166 = t209 + qJ(3,1);
t140 = cos(t166);
t190 = cos(qJ(3,1));
t292 = t140 + t190;
t388 = t292 * pkin(2);
t394 = t399 / 0.2e1 + t388 / 0.2e1;
t165 = t209 + qJ(3,2);
t139 = cos(t165);
t188 = cos(qJ(3,2));
t293 = t139 + t188;
t387 = t293 * pkin(2);
t393 = t398 / 0.2e1 + t387 / 0.2e1;
t164 = t209 + qJ(3,3);
t138 = cos(t164);
t186 = cos(qJ(3,3));
t294 = t138 + t186;
t386 = t294 * pkin(2);
t392 = t397 / 0.2e1 + t386 / 0.2e1;
t199 = xDP(3);
t170 = t199 ^ 2;
t218 = 0.1e1 / pkin(3);
t385 = t170 * t218;
t347 = pkin(5) + qJ(2,1);
t157 = (rSges(3,3) + t347);
t357 = m(3) * t157;
t96 = t357 + m(2) * (rSges(2,3) + qJ(2,1));
t346 = pkin(5) + qJ(2,2);
t156 = (rSges(3,3) + t346);
t358 = m(3) * t156;
t95 = t358 + m(2) * (rSges(2,3) + qJ(2,2));
t345 = pkin(5) + qJ(2,3);
t155 = (rSges(3,3) + t345);
t359 = m(3) * t155;
t94 = t359 + m(2) * (rSges(2,3) + qJ(2,3));
t197 = rSges(2,2) * m(2);
t325 = sin(pkin(7));
t326 = cos(pkin(7));
t368 = m(2) * rSges(2,1);
t376 = pkin(2) * m(3);
t229 = -t325 * t197 + (t368 + t376) * t326;
t383 = 4 * rSges(2,3);
t141 = cos(t167);
t126 = 0.1e1 / t141;
t142 = cos(t168);
t128 = 0.1e1 / t142;
t143 = cos(t169);
t130 = 0.1e1 / t143;
t382 = 0.1e1 / t141 ^ 2;
t381 = 0.1e1 / t142 ^ 2;
t380 = 0.1e1 / t143 ^ 2;
t366 = m(3) * rSges(3,2);
t124 = rSges(3,1) * t366 - Icges(3,4);
t379 = 0.2e1 * t124;
t378 = m(2) / 0.2e1;
t377 = m(3) / 0.2e1;
t108 = sin(t120);
t135 = sin(t167);
t132 = sin(t164);
t180 = sin(qJ(3,3));
t297 = t132 + t180;
t61 = t297 * pkin(2) + pkin(3) * t108 + t135 * t396;
t375 = -t61 / 0.2e1;
t374 = t61 / 0.2e1;
t109 = sin(t121);
t136 = sin(t168);
t133 = sin(t165);
t182 = sin(qJ(3,2));
t296 = t133 + t182;
t62 = t296 * pkin(2) + pkin(3) * t109 + t136 * t396;
t373 = -t62 / 0.2e1;
t372 = t62 / 0.2e1;
t110 = sin(t122);
t137 = sin(t169);
t134 = sin(t166);
t184 = sin(qJ(3,1));
t295 = t134 + t184;
t63 = t295 * pkin(2) + pkin(3) * t110 + t137 * t396;
t371 = -t63 / 0.2e1;
t370 = t63 / 0.2e1;
t213 = rSges(3,2) ^ 2;
t215 = rSges(3,1) ^ 2;
t99 = m(3) * (-t213 + t215) - Icges(3,1) + Icges(3,2);
t369 = t99 / 0.2e1;
t367 = m(3) * rSges(3,1);
t365 = rSges(3,1) * g(3);
t364 = -t199 / 0.2e1;
t363 = t202 / 0.2e1;
t356 = pkin(3) * t141;
t355 = pkin(3) * t142;
t354 = pkin(3) * t143;
t127 = t126 * t382;
t158 = -pkin(6) - t345;
t151 = 0.1e1 / t158;
t321 = t126 * t199;
t287 = t61 * t321;
t200 = xDP(2);
t174 = legFrame(3,3);
t145 = sin(t174);
t148 = cos(t174);
t144 = t326 * pkin(2);
t119 = t144 + pkin(1);
t181 = sin(qJ(1,3));
t187 = cos(qJ(1,3));
t64 = t119 * t181 + t158 * t187;
t67 = t119 * t187 - t158 * t181;
t77 = t145 * t187 + t148 * t181;
t46 = t145 * t67 + t148 * t64 + t77 * t356;
t329 = t46 * t200;
t201 = xDP(1);
t74 = -t145 * t181 + t148 * t187;
t43 = -t145 * t64 + t148 * t67 + t74 * t356;
t332 = t43 * t201;
t235 = t287 + 0.2e1 * t329 + 0.2e1 * t332;
t323 = t126 * t151;
t274 = -t323 / 0.2e1;
t28 = (t332 + t329 + t287 / 0.2e1) * t151;
t52 = (t135 * t321 + t200 * t77 + t201 * t74) * t151;
t344 = t126 * t52;
t49 = pkin(1) * t52;
t10 = (-t127 * t385 - (-t28 + t344 * t392 + (t49 * t126 + t235 * t274) * t141) * t52) * t151;
t101 = -rSges(3,2) * t359 + Icges(3,6);
t275 = m(1) * rSges(1,1) + t400;
t338 = t186 * rSges(3,2);
t341 = t180 * rSges(3,1);
t232 = -(t197 + (t338 + t341) * m(3)) * t325 - (-t368 + (-rSges(3,1) * t186 + rSges(3,2) * t180 - pkin(2)) * m(3)) * t326 + t275;
t303 = t199 * t218;
t288 = m(3) * t303;
t267 = pkin(2) * t288;
t251 = t126 * t267;
t255 = rSges(3,1) * t359 - Icges(3,5);
t304 = t170 / pkin(3) ^ 2;
t262 = t127 * t135 * t304;
t302 = t218 * t126;
t281 = t199 * t302;
t265 = t108 * t281;
t266 = -0.2e1 * t124 * t303;
t198 = m(1) * rSges(1,2);
t270 = t198 - t94;
t283 = -0.2e1 * pkin(1) * t366;
t284 = t367 * t396;
t214 = rSges(2,2) ^ 2;
t216 = rSges(2,1) ^ 2;
t220 = pkin(2) ^ 2;
t228 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (-rSges(2,1) * t197 + Icges(2,4)) * sin(t209) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 + (m(3) * t220 + (-t214 + t216) * m(2) - Icges(2,1) + Icges(2,2)) * cos(t209) / 0.2e1;
t241 = -t135 * t366 + t141 * t367;
t221 = pkin(1) ^ 2;
t206 = 2 * t221;
t291 = t213 + t215;
t258 = t206 + t220 + t291;
t259 = (2 * rSges(2,3) ^ 2) + t206 + t214 + t216;
t31 = t111 * t369 + ((2 * t155 ^ 2) + t258) * t377 + (((t383 + 2 * qJ(2,3)) * qJ(2,3)) + t259) * t378 - t124 * t108 + (t294 * rSges(3,1) - t297 * rSges(3,2)) * t376 + (t229 + t241) * t396 + t228;
t252 = t229 + t400;
t55 = t241 + t252;
t58 = t101 * t141 - t255 * t135;
t208 = 3 * pkin(7);
t210 = 2 * qJ(3,3);
t217 = pkin(3) ^ 2;
t248 = -0.4e1 * pkin(3) * t144;
t276 = -0.3e1 * t217 - 0.2e1 * t220 - (4 * t221);
t298 = -0.2e1 * pkin(2) * pkin(3);
t7 = (t119 + t356) * t382 * t385 * t323 + (t158 * t108 * t382 * t364 + (-(-t217 * cos((3 * t167)) + (-0.4e1 * t158 ^ 2 + t276) * t141 + t248 + (cos((t210 + t208)) + cos((pkin(7) + t210))) * t298 + (-cos((-pkin(7) + qJ(3,3))) - cos((t208 + qJ(3,3)))) * t220) * t52 / 0.4e1 + t158 * t265 * t395 + (-t386 - t397) * (-t49 - (-t332 / 0.2e1 - t329 / 0.2e1 - t287 / 0.4e1) * t151) - (t141 * t396 + t392) * t28) * t126) * t151 * t52;
t87 = -g(1) * t145 + g(2) * t148;
t90 = g(1) * t148 + g(2) * t145;
t1 = t31 * t10 + t55 * t7 + t58 * t262 + t52 * t99 * t265 - t111 * t266 * t344 + (-t255 * t281 - t52 * t283) * t303 - (t101 * t281 - t52 * t284) * t135 * t281 - 0.2e1 * t52 * (-t28 * t94 + (-t338 / 0.2e1 - t341 / 0.2e1) * t251) + (-t232 * t87 + t270 * t90) * t187 + (t232 * t90 + t270 * t87) * t181 - (-rSges(3,1) * t132 - rSges(3,2) * t138) * t52 * t251;
t353 = t151 * t1;
t271 = -0.2e1 * t288;
t118 = rSges(3,2) * t271;
t257 = rSges(3,1) * t271;
t4 = -t55 * t10 + t52 * (t126 * t135 * t257 - t52 * t94 + t118) + (-t181 * t90 + t187 * t87 - t7) * t202;
t352 = t151 * t4;
t159 = -pkin(6) - t346;
t152 = 0.1e1 / t159;
t102 = -rSges(3,2) * t358 + Icges(3,6);
t129 = t128 * t381;
t318 = t128 * t199;
t286 = t62 * t318;
t175 = legFrame(2,3);
t146 = sin(t175);
t149 = cos(t175);
t183 = sin(qJ(1,2));
t189 = cos(qJ(1,2));
t65 = t119 * t183 + t159 * t189;
t68 = t119 * t189 - t159 * t183;
t78 = t146 * t189 + t149 * t183;
t47 = t146 * t68 + t149 * t65 + t78 * t355;
t328 = t47 * t200;
t75 = -t146 * t183 + t149 * t189;
t44 = -t146 * t65 + t149 * t68 + t75 * t355;
t331 = t44 * t201;
t234 = t286 + 0.2e1 * t328 + 0.2e1 * t331;
t320 = t128 * t152;
t273 = -t320 / 0.2e1;
t29 = (t331 + t328 + t286 / 0.2e1) * t152;
t53 = (t136 * t318 + t200 * t78 + t201 * t75) * t152;
t343 = t128 * t53;
t50 = pkin(1) * t53;
t11 = (-t129 * t385 - (-t29 + t343 * t393 + (t50 * t128 + t234 * t273) * t142) * t53) * t152;
t337 = t188 * rSges(3,2);
t340 = t182 * rSges(3,1);
t231 = -(t197 + (t337 + t340) * m(3)) * t325 - (-t368 + (-rSges(3,1) * t188 + rSges(3,2) * t182 - pkin(2)) * m(3)) * t326 + t275;
t250 = t128 * t267;
t254 = rSges(3,1) * t358 - Icges(3,5);
t261 = t129 * t136 * t304;
t301 = t218 * t128;
t279 = t199 * t301;
t264 = t109 * t279;
t269 = t198 - t95;
t240 = -t136 * t366 + t142 * t367;
t32 = t112 * t369 + ((2 * t156 ^ 2) + t258) * t377 + (((t383 + 2 * qJ(2,2)) * qJ(2,2)) + t259) * t378 - t124 * t109 + (t293 * rSges(3,1) - t296 * rSges(3,2)) * t376 + (t229 + t240) * t396 + t228;
t56 = t240 + t252;
t59 = t102 * t142 - t254 * t136;
t211 = 2 * qJ(3,2);
t8 = (t119 + t355) * t381 * t385 * t320 + (t159 * t109 * t381 * t364 + (-(-t217 * cos((3 * t168)) + (-0.4e1 * t159 ^ 2 + t276) * t142 + t248 + (cos((t208 + t211)) + cos((pkin(7) + t211))) * t298 + (-cos((-pkin(7) + qJ(3,2))) - cos((t208 + qJ(3,2)))) * t220) * t53 / 0.4e1 + t159 * t264 * t395 + (-t387 - t398) * (-t50 - (-t331 / 0.2e1 - t328 / 0.2e1 - t286 / 0.4e1) * t152) - (t142 * t396 + t393) * t29) * t128) * t152 * t53;
t88 = -g(1) * t146 + g(2) * t149;
t91 = g(1) * t149 + g(2) * t146;
t2 = t32 * t11 + t56 * t8 + t59 * t261 + t53 * t99 * t264 - t112 * t266 * t343 + (-t254 * t279 - t53 * t283) * t303 - (t102 * t279 - t53 * t284) * t136 * t279 - 0.2e1 * t53 * (-t29 * t95 + (-t337 / 0.2e1 - t340 / 0.2e1) * t250) + (-t231 * t88 + t269 * t91) * t189 + (t231 * t91 + t269 * t88) * t183 - (-rSges(3,1) * t133 - rSges(3,2) * t139) * t53 * t250;
t351 = t152 * t2;
t5 = -t56 * t11 + t53 * (t128 * t136 * t257 - t53 * t95 + t118) + (-t183 * t91 + t189 * t88 - t8) * t202;
t350 = t152 * t5;
t160 = -pkin(6) - t347;
t153 = 0.1e1 / t160;
t103 = -rSges(3,2) * t357 + Icges(3,6);
t131 = t130 * t380;
t315 = t130 * t199;
t285 = t63 * t315;
t176 = legFrame(1,3);
t147 = sin(t176);
t150 = cos(t176);
t185 = sin(qJ(1,1));
t191 = cos(qJ(1,1));
t66 = t119 * t185 + t160 * t191;
t69 = t119 * t191 - t160 * t185;
t79 = t147 * t191 + t150 * t185;
t48 = t147 * t69 + t150 * t66 + t79 * t354;
t327 = t48 * t200;
t76 = -t147 * t185 + t150 * t191;
t45 = -t147 * t66 + t150 * t69 + t76 * t354;
t330 = t45 * t201;
t233 = t285 + 0.2e1 * t327 + 0.2e1 * t330;
t317 = t130 * t153;
t272 = -t317 / 0.2e1;
t30 = (t330 + t327 + t285 / 0.2e1) * t153;
t54 = (t137 * t315 + t200 * t79 + t201 * t76) * t153;
t342 = t130 * t54;
t51 = pkin(1) * t54;
t12 = (-t131 * t385 - (-t30 + t342 * t394 + (t51 * t130 + t233 * t272) * t143) * t54) * t153;
t336 = t190 * rSges(3,2);
t339 = t184 * rSges(3,1);
t230 = -(t197 + (t336 + t339) * m(3)) * t325 - (-t368 + (-rSges(3,1) * t190 + rSges(3,2) * t184 - pkin(2)) * m(3)) * t326 + t275;
t249 = t130 * t267;
t253 = rSges(3,1) * t357 - Icges(3,5);
t260 = t131 * t137 * t304;
t300 = t218 * t130;
t277 = t199 * t300;
t263 = t110 * t277;
t268 = t198 - t96;
t239 = -t137 * t366 + t143 * t367;
t33 = t113 * t369 + ((2 * t157 ^ 2) + t258) * t377 + (((t383 + 2 * qJ(2,1)) * qJ(2,1)) + t259) * t378 - t124 * t110 + (t292 * rSges(3,1) - t295 * rSges(3,2)) * t376 + (t229 + t239) * t396 + t228;
t57 = t239 + t252;
t60 = t103 * t143 - t253 * t137;
t89 = -g(1) * t147 + g(2) * t150;
t212 = 2 * qJ(3,1);
t9 = (t119 + t354) * t380 * t385 * t317 + (t160 * t110 * t380 * t364 + (-(-t217 * cos((3 * t169)) + (-0.4e1 * t160 ^ 2 + t276) * t143 + t248 + (cos((t208 + t212)) + cos((pkin(7) + t212))) * t298 + (-cos((-pkin(7) + qJ(3,1))) - cos((t208 + qJ(3,1)))) * t220) * t54 / 0.4e1 + t160 * t263 * t395 + (-t388 - t399) * (-t51 - (-t330 / 0.2e1 - t327 / 0.2e1 - t285 / 0.4e1) * t153) - (t143 * t396 + t394) * t30) * t130) * t153 * t54;
t92 = g(1) * t150 + g(2) * t147;
t3 = t33 * t12 + t57 * t9 + t60 * t260 + t54 * t99 * t263 - t113 * t266 * t342 + (-t253 * t277 - t54 * t283) * t303 - (t103 * t277 - t54 * t284) * t137 * t277 - 0.2e1 * t54 * (-t30 * t96 + (-t336 / 0.2e1 - t339 / 0.2e1) * t249) + (-t230 * t89 + t268 * t92) * t191 + (t230 * t92 + t268 * t89) * t185 - (-rSges(3,1) * t134 - rSges(3,2) * t140) * t54 * t249;
t349 = t153 * t3;
t6 = -t57 * t12 + t54 * (t130 * t137 * t257 - t54 * t96 + t118) + (-t185 * t92 + t191 * t89 - t9) * t202;
t348 = t153 * t6;
t335 = t218 * t58;
t334 = t218 * t59;
t333 = t218 * t60;
t107 = t291 * m(3) + Icges(3,3);
t324 = t107 * t218;
t177 = xDDP(3);
t322 = t126 * t177;
t319 = t128 * t177;
t316 = t130 * t177;
t314 = t135 * t151;
t313 = t136 * t152;
t312 = t137 * t153;
t178 = xDDP(2);
t311 = t151 * t178;
t179 = xDDP(1);
t310 = t151 * t179;
t309 = t152 * t178;
t308 = t152 * t179;
t307 = t153 * t178;
t306 = t153 * t179;
t290 = 0.2e1 * m(3);
t282 = t151 * t322;
t280 = t152 * t319;
t278 = t153 * t316;
t247 = t181 * t87 + t187 * t90;
t246 = t183 * t88 + t189 * t91;
t245 = t185 * t89 + t191 * t92;
t196 = rSges(3,2) * g(3);
t42 = (-t137 * t57 + t63 * t363) * t317;
t41 = (-t136 * t56 + t62 * t363) * t320;
t40 = (-t135 * t55 + t61 * t363) * t323;
t39 = (t202 * t48 - t57 * t79) * t153;
t38 = (t202 * t47 - t56 * t78) * t152;
t37 = (t202 * t46 - t55 * t77) * t151;
t36 = (t202 * t45 - t57 * t76) * t153;
t35 = (t202 * t44 - t56 * t75) * t152;
t34 = (t202 * t43 - t55 * t74) * t151;
t21 = (t333 - (t137 * t33 + t57 * t371) * t153) * t130;
t20 = (t334 - (t136 * t32 + t56 * t373) * t152) * t128;
t19 = (t335 - (t135 * t31 + t55 * t375) * t151) * t126;
t18 = (t33 * t79 - t48 * t57) * t153;
t17 = (t32 * t78 - t47 * t56) * t152;
t16 = (t31 * t77 - t46 * t55) * t151;
t15 = (t33 * t76 - t45 * t57) * t153;
t14 = (t32 * t75 - t44 * t56) * t152;
t13 = (t31 * t74 - t43 * t55) * t151;
t22 = [-(-t15 * t76 - t36 * t45) * t306 - (-t15 * t79 - t36 * t48) * t307 - (-t15 * t137 + t76 * t333 - t36 * t370) * t278 - t76 * t349 - t45 * t348 - (-t14 * t75 - t35 * t44) * t308 - (-t14 * t78 - t35 * t47) * t309 - (-t14 * t136 + t75 * t334 - t35 * t372) * t280 - t75 * t351 - t44 * t350 - (-t13 * t74 - t34 * t43) * t310 - (-t13 * t77 - t34 * t46) * t311 - (-t13 * t135 + t74 * t335 - t34 * t374) * t282 - t74 * t353 - t43 * t352 + (-g(1) + t179) * m(4); -(-t18 * t76 - t39 * t45) * t306 - (-t18 * t79 - t39 * t48) * t307 - (-t18 * t137 + t79 * t333 - t39 * t370) * t278 - t79 * t349 - t48 * t348 - (-t17 * t75 - t38 * t44) * t308 - (-t17 * t78 - t38 * t47) * t309 - (-t17 * t136 + t78 * t334 - t38 * t372) * t280 - t78 * t351 - t47 * t350 - (-t16 * t74 - t37 * t43) * t310 - (-t16 * t77 - t37 * t46) * t311 - (-t16 * t135 + t77 * t335 - t37 * t374) * t282 - t77 * t353 - t46 * t352 + (-g(2) + t178) * m(4); -(t21 * t76 - t42 * t45) * t306 - (t21 * t79 - t42 * t48) * t307 + (-t21 * t312 - t42 * t153 * t371 + (-t60 * t312 + t324) * t300) * t316 - t130 * t3 * t312 + t63 * t6 * t272 + (t60 * t12 + t107 * t260 - t54 * ((rSges(3,1) * t137 + rSges(3,2) * t143) * (t233 * t153 - t51) * t290 - (t99 * t110 + t113 * t379 + (t295 * rSges(3,1) + t292 * rSges(3,2)) * t376) * t54) / 0.2e1 + m(3) * ((t245 * rSges(3,2) - t365) * t143 + t137 * (t245 * rSges(3,1) + t196))) * t300 - (t20 * t75 - t41 * t44) * t308 - (t20 * t78 - t41 * t47) * t309 + (-t20 * t313 - t41 * t152 * t373 + (-t59 * t313 + t324) * t301) * t319 - t128 * t2 * t313 + t62 * t5 * t273 + (t59 * t11 + t107 * t261 - t53 * ((rSges(3,1) * t136 + rSges(3,2) * t142) * (t234 * t152 - t50) * t290 - (t99 * t109 + t112 * t379 + (t296 * rSges(3,1) + t293 * rSges(3,2)) * t376) * t53) / 0.2e1 + m(3) * ((t246 * rSges(3,2) - t365) * t142 + t136 * (t246 * rSges(3,1) + t196))) * t301 - (t19 * t74 - t40 * t43) * t310 - (t19 * t77 - t40 * t46) * t311 + (-t19 * t314 - t40 * t151 * t375 + (-t58 * t314 + t324) * t302) * t322 - t126 * t1 * t314 + t61 * t4 * t274 + (t58 * t10 + t107 * t262 - t52 * ((rSges(3,1) * t135 + rSges(3,2) * t141) * (t235 * t151 - t49) * t290 - (t99 * t108 + t111 * t379 + (t297 * rSges(3,1) + t294 * rSges(3,2)) * t376) * t52) / 0.2e1 + m(3) * ((t247 * rSges(3,2) - t365) * t141 + t135 * (t247 * rSges(3,1) + t196))) * t302 + (-g(3) + t177) * m(4);];
tauX  = t22;
