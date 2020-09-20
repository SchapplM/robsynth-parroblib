% Calculate minimal parameter regressor of inverse dynamics forces for
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x15]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR9V1G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:47:17
% EndTime: 2020-08-06 18:47:34
% DurationCPUTime: 16.69s
% Computational Cost: add. (110969->741), mult. (110922->1150), div. (9828->14), fcn. (53133->107), ass. (0->519)
t270 = (pkin(7) + qJ(3,1));
t208 = 2 * t270;
t195 = sin(t208);
t306 = 2 * pkin(7);
t267 = t306 + qJ(3,1);
t226 = sin(t267);
t290 = sin(qJ(3,1));
t459 = t226 + t290;
t229 = sin(t270);
t598 = 0.2e1 * t229;
t111 = pkin(1) * t598 + t459 * pkin(2) + pkin(3) * t195;
t235 = cos(t270);
t221 = 0.1e1 / t235;
t283 = xDDP(3);
t417 = t111 * t221 * t283;
t285 = xDDP(1);
t276 = cos(pkin(7));
t205 = t276 * pkin(2) + pkin(1);
t282 = pkin(5) + qJ(2,1);
t260 = -pkin(6) - t282;
t291 = sin(qJ(1,1));
t297 = cos(qJ(1,1));
t123 = t291 * t205 + t297 * t260;
t126 = t205 * t297 - t291 * t260;
t279 = legFrame(1,3);
t243 = sin(t279);
t246 = cos(t279);
t360 = t243 * t291 - t246 * t297;
t557 = pkin(3) * t235;
t84 = -t123 * t243 + t126 * t246 - t360 * t557;
t509 = t285 * t84;
t284 = xDDP(2);
t135 = t243 * t297 + t246 * t291;
t87 = t123 * t246 + t126 * t243 + t135 * t557;
t512 = t284 * t87;
t629 = t417 / 0.2e1 + t509 + t512;
t269 = pkin(7) + qJ(3,2);
t207 = 2 * t269;
t194 = sin(t207);
t266 = t306 + qJ(3,2);
t225 = sin(t266);
t288 = sin(qJ(3,2));
t460 = t225 + t288;
t228 = sin(t269);
t599 = 0.2e1 * t228;
t110 = pkin(1) * t599 + t460 * pkin(2) + pkin(3) * t194;
t234 = cos(t269);
t218 = 0.1e1 / t234;
t418 = t110 * t218 * t283;
t281 = pkin(5) + qJ(2,2);
t259 = -pkin(6) - t281;
t289 = sin(qJ(1,2));
t295 = cos(qJ(1,2));
t122 = t289 * t205 + t295 * t259;
t125 = t205 * t295 - t289 * t259;
t278 = legFrame(2,3);
t242 = sin(t278);
t245 = cos(t278);
t361 = t242 * t289 - t245 * t295;
t558 = pkin(3) * t234;
t83 = -t122 * t242 + t125 * t245 - t361 * t558;
t510 = t285 * t83;
t134 = t242 * t295 + t245 * t289;
t86 = t122 * t245 + t125 * t242 + t134 * t558;
t513 = t284 * t86;
t628 = t418 / 0.2e1 + t510 + t513;
t268 = pkin(7) + qJ(3,3);
t206 = 2 * t268;
t193 = sin(t206);
t265 = t306 + qJ(3,3);
t224 = sin(t265);
t286 = sin(qJ(3,3));
t461 = t224 + t286;
t227 = sin(t268);
t600 = 0.2e1 * t227;
t109 = pkin(1) * t600 + t461 * pkin(2) + pkin(3) * t193;
t233 = cos(t268);
t215 = 0.1e1 / t233;
t419 = t109 * t215 * t283;
t280 = pkin(5) + qJ(2,3);
t258 = -pkin(6) - t280;
t287 = sin(qJ(1,3));
t293 = cos(qJ(1,3));
t121 = t287 * t205 + t293 * t258;
t124 = t205 * t293 - t287 * t258;
t277 = legFrame(3,3);
t241 = sin(t277);
t244 = cos(t277);
t362 = t241 * t287 - t244 * t293;
t559 = pkin(3) * t233;
t82 = -t121 * t241 + t124 * t244 - t362 * t559;
t511 = t285 * t82;
t133 = t241 * t293 + t244 * t287;
t85 = t121 * t244 + t124 * t241 + t133 * t559;
t514 = t284 * t85;
t627 = t419 / 0.2e1 + t511 + t514;
t216 = 0.1e1 / t233 ^ 2;
t248 = 0.1e1 / t258;
t321 = t258 ^ 2;
t249 = 0.1e1 / t321;
t298 = xDP(3);
t312 = 0.1e1 / pkin(3);
t471 = t298 * t312;
t299 = xDP(2);
t300 = xDP(1);
t91 = ((t287 * t299 + t293 * t300) * t244 + t241 * (-t287 * t300 + t293 * t299)) * t233 + t227 * t298;
t508 = t298 * (-t91 * t193 / 0.2e1 + (t205 + t559) * t471) * t215;
t425 = t248 * t508;
t505 = t362 * t248;
t115 = t285 * t505;
t502 = t133 * t248;
t118 = t284 * t502;
t487 = t215 * t248;
t416 = t227 * t487;
t127 = t283 * t416;
t217 = t215 * t216;
t274 = t298 ^ 2;
t481 = t274 * t312;
t139 = t217 * t248 * t481;
t250 = t248 * t249;
t525 = t216 * t91;
t405 = t525 / 0.4e1;
t486 = t216 * t249;
t439 = t91 * t486;
t196 = cos(t206);
t292 = cos(qJ(3,3));
t440 = t91 * t487;
t392 = pkin(3) * t440;
t415 = t233 * t487;
t230 = cos(t265);
t562 = pkin(2) * t230;
t582 = t91 * pkin(1);
t236 = t277 + qJ(1,3);
t201 = qJ(3,3) + t236;
t181 = pkin(7) + t201;
t163 = sin(t181);
t202 = -qJ(3,3) + t236;
t182 = -pkin(7) + t202;
t164 = sin(t182);
t169 = cos(t181);
t170 = cos(t182);
t308 = 2 * qJ(3,3);
t175 = t306 + t308 + t236;
t176 = t306 + t201;
t307 = -2 * pkin(7);
t177 = t307 + t202;
t178 = t307 - (2 * qJ(3,3)) + t236;
t607 = 0.2e1 * pkin(1);
t239 = t299 * t607;
t257 = pkin(1) * t300;
t623 = 0.4e1 * pkin(1);
t448 = t298 * t623;
t592 = -0.2e1 * t300;
t620 = 0.2e1 * t298;
t67 = t227 * t448 + (t163 + t164) * (t258 * t592 + t239) + 0.2e1 * (t169 + t170) * (t258 * t299 + t257) + (t193 * t620 + (cos(t178) + cos(t175) + 0.2e1 * cos(t236)) * t300 + (sin(t178) + sin(t175) + 0.2e1 * sin(t236)) * t299) * pkin(3) + (t461 * t620 + (cos(t177) + cos(t176) + cos(t202) + cos(t201)) * t300 + (sin(t177) + sin(t176) + sin(t202) + sin(t201)) * t299) * pkin(2);
t585 = t67 / 0.2e1;
t55 = t196 * t392 - (-0.2e1 * t582 + t585) * t415 + (t292 * pkin(2) + pkin(3) + t562) * t440;
t46 = t115 - t118 - t127 - t55 * t439 / 0.2e1 + t250 * t67 * t405 - t139;
t43 = t46 * pkin(1);
t480 = t274 / pkin(3) ^ 2;
t446 = -0.2e1 * t480;
t616 = (t230 + t292) * pkin(2);
t65 = (t607 + ((t196 + 0.1e1) * pkin(3) + t616) * t215) * t91 * t248;
t534 = t65 * t67;
t305 = 3 * pkin(7);
t311 = pkin(3) ^ 2;
t314 = pkin(2) ^ 2;
t377 = pkin(2) * t392;
t396 = 0.2e1 * pkin(3) * t471;
t400 = -0.4e1 * pkin(1) ^ 2 - 0.3e1 * t311 - 0.2e1 * t314;
t603 = 0.4e1 * (t582 - t67 / 0.8e1) * t487;
t455 = pkin(3) * t603;
t597 = 0.4e1 * t276;
t28 = -t258 * t215 * t193 * t396 + t196 * t455 - ((-0.4e1 * t321 + t400) * t91 + pkin(1) * t67) * t415 + t377 * t597 + t455 + t603 * t616 + 0.2e1 * (cos((t308 + t305)) + cos((t308 + pkin(7)))) * t377 + (t311 * cos((3 * t268)) + (cos((qJ(3,3) - pkin(7))) + cos((t305 + qJ(3,3)))) * t314) * t440;
t544 = t28 * t91;
t626 = t216 * (t249 * (t544 / 0.2e1 - t534 / 0.4e1) + t280 * t446 + 0.2e1 * t425) + (t419 + 0.2e1 * t511 + 0.2e1 * t514) * t248 + 0.4e1 * t43;
t219 = 0.1e1 / t234 ^ 2;
t251 = 0.1e1 / t259;
t323 = t259 ^ 2;
t252 = 0.1e1 / t323;
t92 = ((t289 * t299 + t295 * t300) * t245 + t242 * (-t289 * t300 + t295 * t299)) * t234 + t228 * t298;
t507 = t298 * (-t92 * t194 / 0.2e1 + (t205 + t558) * t471) * t218;
t424 = t251 * t507;
t504 = t361 * t251;
t116 = t285 * t504;
t501 = t134 * t251;
t119 = t284 * t501;
t485 = t218 * t251;
t414 = t228 * t485;
t128 = t283 * t414;
t220 = t218 * t219;
t140 = t220 * t251 * t481;
t253 = t251 * t252;
t523 = t219 * t92;
t404 = t523 / 0.4e1;
t484 = t219 * t252;
t434 = t92 * t484;
t197 = cos(t207);
t294 = cos(qJ(3,2));
t435 = t92 * t485;
t391 = pkin(3) * t435;
t413 = t234 * t485;
t231 = cos(t266);
t561 = pkin(2) * t231;
t581 = t92 * pkin(1);
t237 = t278 + qJ(1,2);
t203 = qJ(3,2) + t237;
t185 = pkin(7) + t203;
t165 = sin(t185);
t204 = -qJ(3,2) + t237;
t186 = -pkin(7) + t204;
t166 = sin(t186);
t171 = cos(t185);
t172 = cos(t186);
t309 = 2 * qJ(3,2);
t402 = t306 + t237;
t183 = t309 + t402;
t401 = t307 + t237;
t184 = -(2 * qJ(3,2)) + t401;
t187 = qJ(3,2) + t402;
t188 = -qJ(3,2) + t401;
t68 = t228 * t448 + (t165 + t166) * (t259 * t592 + t239) + 0.2e1 * (t171 + t172) * (t259 * t299 + t257) + (t194 * t620 + (cos(t184) + cos(t183) + 0.2e1 * cos(t237)) * t300 + (sin(t184) + sin(t183) + 0.2e1 * sin(t237)) * t299) * pkin(3) + (t460 * t620 + (cos(t188) + cos(t187) + cos(t204) + cos(t203)) * t300 + (sin(t188) + sin(t187) + sin(t204) + sin(t203)) * t299) * pkin(2);
t584 = t68 / 0.2e1;
t56 = t197 * t391 - (-0.2e1 * t581 + t584) * t413 + (t294 * pkin(2) + pkin(3) + t561) * t435;
t47 = t116 - t119 - t128 - t56 * t434 / 0.2e1 + t253 * t68 * t404 - t140;
t44 = t47 * pkin(1);
t615 = (t231 + t294) * pkin(2);
t66 = (t607 + ((t197 + 0.1e1) * pkin(3) + t615) * t218) * t92 * t251;
t533 = t66 * t68;
t376 = pkin(2) * t391;
t602 = 0.4e1 * (t581 - t68 / 0.8e1) * t485;
t454 = pkin(3) * t602;
t29 = -t259 * t218 * t194 * t396 + t197 * t454 - ((-0.4e1 * t323 + t400) * t92 + pkin(1) * t68) * t413 + t376 * t597 + t454 + t602 * t615 + 0.2e1 * (cos((t309 + t305)) + cos((t309 + pkin(7)))) * t376 + (t311 * cos((3 * t269)) + (cos((qJ(3,2) - pkin(7))) + cos((t305 + qJ(3,2)))) * t314) * t435;
t543 = t29 * t92;
t625 = t219 * (t252 * (t543 / 0.2e1 - t533 / 0.4e1) + t281 * t446 + 0.2e1 * t424) + (t418 + 0.2e1 * t510 + 0.2e1 * t513) * t251 + 0.4e1 * t44;
t222 = 0.1e1 / t235 ^ 2;
t254 = 0.1e1 / t260;
t325 = t260 ^ 2;
t255 = 0.1e1 / t325;
t93 = ((t291 * t299 + t297 * t300) * t246 + t243 * (-t291 * t300 + t297 * t299)) * t235 + t229 * t298;
t506 = t298 * (-t93 * t195 / 0.2e1 + (t205 + t557) * t471) * t221;
t423 = t254 * t506;
t503 = t360 * t254;
t117 = t285 * t503;
t500 = t135 * t254;
t120 = t284 * t500;
t483 = t221 * t254;
t412 = t229 * t483;
t129 = t283 * t412;
t223 = t221 * t222;
t141 = t223 * t254 * t481;
t256 = t254 * t255;
t521 = t222 * t93;
t403 = t521 / 0.4e1;
t482 = t222 * t255;
t429 = t93 * t482;
t198 = cos(t208);
t296 = cos(qJ(3,1));
t430 = t93 * t483;
t390 = pkin(3) * t430;
t411 = t235 * t483;
t232 = cos(t267);
t560 = pkin(2) * t232;
t580 = t93 * pkin(1);
t238 = t279 + qJ(1,1);
t199 = qJ(3,1) + t238;
t189 = pkin(7) + t199;
t167 = sin(t189);
t200 = -qJ(3,1) + t238;
t190 = -pkin(7) + t200;
t168 = sin(t190);
t173 = cos(t189);
t174 = cos(t190);
t179 = t306 + t199;
t180 = t307 + t200;
t310 = 2 * qJ(3,1);
t191 = t306 + t310 + t238;
t192 = t307 - (2 * qJ(3,1)) + t238;
t69 = t229 * t448 + (t167 + t168) * (t260 * t592 + t239) + 0.2e1 * (t173 + t174) * (t260 * t299 + t257) + (t195 * t620 + (cos(t192) + cos(t191) + 0.2e1 * cos(t238)) * t300 + (sin(t192) + sin(t191) + 0.2e1 * sin(t238)) * t299) * pkin(3) + (t459 * t620 + (cos(t180) + cos(t179) + cos(t200) + cos(t199)) * t300 + (sin(t180) + sin(t179) + sin(t200) + sin(t199)) * t299) * pkin(2);
t583 = t69 / 0.2e1;
t57 = t198 * t390 - (-0.2e1 * t580 + t583) * t411 + (pkin(2) * t296 + pkin(3) + t560) * t430;
t48 = t117 - t120 - t129 - t57 * t429 / 0.2e1 + t256 * t69 * t403 - t141;
t45 = t48 * pkin(1);
t614 = (t232 + t296) * pkin(2);
t64 = (t607 + ((t198 + 0.1e1) * pkin(3) + t614) * t221) * t93 * t254;
t535 = t64 * t69;
t375 = pkin(2) * t390;
t601 = 0.4e1 * (t580 - t69 / 0.8e1) * t483;
t453 = pkin(3) * t601;
t30 = -t260 * t221 * t195 * t396 + t198 * t453 - ((-0.4e1 * t325 + t400) * t93 + pkin(1) * t69) * t411 + t375 * t597 + t453 + t601 * t614 + 0.2e1 * (cos((t305 + t310)) + cos((pkin(7) + t310))) * t375 + (t311 * cos((3 * t270)) + (cos((qJ(3,1) - pkin(7))) + cos((qJ(3,1) + t305))) * t314) * t430;
t542 = t30 * t93;
t624 = t222 * (t255 * (t542 / 0.2e1 - t535 / 0.4e1) + t282 * t446 + 0.2e1 * t423) + (t417 + 0.2e1 * t509 + 0.2e1 * t512) * t254 + 0.4e1 * t45;
t275 = sin(pkin(7));
t622 = 0.2e1 * t275;
t621 = -0.2e1 * t276;
t427 = t222 * t506;
t619 = t30 * t255 * t403 - t482 * t535 / 0.8e1 + (t427 + t629) * t254;
t432 = t219 * t507;
t618 = t29 * t252 * t404 - t484 * t533 / 0.8e1 + (t432 + t628) * t251;
t437 = t216 * t508;
t617 = t28 * t249 * t405 - t486 * t534 / 0.8e1 + (t437 + t627) * t248;
t545 = t246 * g(2);
t552 = t243 * g(1);
t144 = -t545 + t552;
t546 = t246 * g(1);
t551 = t243 * g(2);
t147 = t546 + t551;
t108 = -t144 * t291 + t147 * t297;
t105 = t144 * t297 + t147 * t291;
t547 = t245 * g(2);
t554 = t242 * g(1);
t143 = -t547 + t554;
t548 = t245 * g(1);
t553 = t242 * g(2);
t146 = t548 + t553;
t107 = -t143 * t289 + t146 * t295;
t104 = t143 * t295 + t146 * t289;
t549 = t244 * g(2);
t556 = t241 * g(1);
t142 = -t549 + t556;
t550 = t244 * g(1);
t555 = t241 * g(2);
t145 = t550 + t555;
t106 = -t142 * t287 + t145 * t293;
t103 = t142 * t293 + t145 * t287;
t399 = t67 / 0.16e2;
t19 = (t550 / 0.2e1 + t555 / 0.2e1) * t287 + t43 + (-t549 / 0.2e1 + t556 / 0.2e1) * t293 + (t544 / 0.8e1 - t65 * t399) * t486 - (-t511 / 0.2e1 - t514 / 0.2e1 - t419 / 0.4e1 - t437 / 0.2e1) * t248;
t606 = 0.2e1 * t19;
t398 = t68 / 0.16e2;
t20 = (t548 / 0.2e1 + t553 / 0.2e1) * t289 + t44 + (-t547 / 0.2e1 + t554 / 0.2e1) * t295 + (t543 / 0.8e1 - t66 * t398) * t484 - (-t510 / 0.2e1 - t513 / 0.2e1 - t418 / 0.4e1 - t432 / 0.2e1) * t251;
t605 = 0.2e1 * t20;
t397 = t69 / 0.16e2;
t21 = (t546 / 0.2e1 + t551 / 0.2e1) * t291 + t45 + (-t545 / 0.2e1 + t552 / 0.2e1) * t297 + (t542 / 0.8e1 - t64 * t397) * t482 - (-t509 / 0.2e1 - t512 / 0.2e1 - t417 / 0.4e1 - t427 / 0.2e1) * t254;
t604 = 0.2e1 * t21;
t596 = -0.2e1 * t280;
t595 = -0.2e1 * t281;
t594 = -0.2e1 * t282;
t590 = -g(1) / 0.2e1;
t589 = g(1) / 0.2e1;
t588 = -g(2) / 0.2e1;
t587 = g(2) / 0.2e1;
t586 = pkin(1) * g(2);
t579 = t109 / 0.2e1;
t578 = t110 / 0.2e1;
t577 = t111 / 0.2e1;
t576 = t163 / 0.2e1;
t575 = t165 / 0.2e1;
t574 = t167 / 0.2e1;
t573 = t170 / 0.2e1;
t572 = t172 / 0.2e1;
t571 = t174 / 0.2e1;
t570 = t227 / 0.2e1;
t569 = t228 / 0.2e1;
t568 = t229 / 0.2e1;
t567 = t233 / 0.2e1;
t566 = t234 / 0.2e1;
t565 = t235 / 0.2e1;
t264 = t276 ^ 2;
t564 = -0.1e1 + 0.2e1 * t264;
t563 = 0.4e1 * t264 - 0.2e1;
t541 = t46 * t82;
t540 = t46 * t85;
t539 = t47 * t83;
t538 = t47 * t86;
t537 = t48 * t84;
t536 = t48 * t87;
t532 = qJ(2,1) * t48;
t531 = qJ(2,2) * t47;
t530 = qJ(2,3) * t46;
t529 = t109 * t46;
t528 = t110 * t47;
t527 = t111 * t48;
t526 = t215 * t46;
t524 = t218 * t47;
t522 = t221 * t48;
t88 = t91 ^ 2;
t520 = t250 * t88;
t89 = t92 ^ 2;
t519 = t253 * t89;
t90 = t93 ^ 2;
t518 = t256 * t90;
t271 = t292 ^ 2;
t517 = t271 * t46;
t272 = t294 ^ 2;
t516 = t272 * t47;
t273 = t296 ^ 2;
t515 = t273 * t48;
t479 = t275 * t276;
t478 = t286 * t275;
t477 = t286 * t292;
t476 = t290 * t275;
t475 = t290 * t296;
t474 = t292 * t276;
t473 = t294 * t288;
t472 = t296 * t276;
t470 = t312 * t215;
t469 = t312 * t218;
t468 = t312 * t221;
t467 = g(2) * t576 + t169 * t589;
t466 = g(2) * t575 + t171 * t589;
t465 = g(2) * t574 + t173 * t589;
t464 = g(1) * t576 + t169 * t588;
t463 = g(1) * t575 + t171 * t588;
t462 = g(1) * t574 + t173 * t588;
t452 = t264 - 0.1e1 / 0.2e1;
t451 = t271 - 0.1e1 / 0.2e1;
t450 = t272 - 0.1e1 / 0.2e1;
t449 = t273 - 0.1e1 / 0.2e1;
t445 = -0.8e1 * t479;
t444 = 0.4e1 * t479;
t438 = t216 * t520;
t436 = t217 * t249 * t88;
t433 = t219 * t519;
t431 = t220 * t252 * t89;
t428 = t222 * t518;
t426 = t223 * t255 * t90;
t422 = t46 * t477;
t421 = t48 * t475;
t420 = t47 * t473;
t410 = t251 * t471;
t409 = t254 * t471;
t408 = t275 * t480;
t407 = t276 * t480;
t406 = t248 * t471;
t395 = pkin(2) * t88 * t486;
t394 = pkin(2) * t89 * t484;
t393 = pkin(2) * t90 * t482;
t386 = t474 * t478;
t385 = t472 * t476;
t384 = t473 * t479;
t383 = -t109 * t487 / 0.2e1;
t381 = -t110 * t485 / 0.2e1;
t379 = -t111 * t483 / 0.2e1;
t374 = t406 * t525;
t373 = t410 * t523;
t372 = t409 * t521;
t371 = t439 * t585;
t370 = t434 * t584;
t369 = t429 * t583;
t359 = pkin(2) * t374;
t358 = pkin(2) * t373;
t357 = pkin(2) * t372;
t356 = g(2) * t573 + t164 * t590;
t355 = g(2) * t572 + t166 * t590;
t354 = g(2) * t571 + t168 * t590;
t353 = g(1) * t573 + t164 * t587;
t352 = g(1) * t572 + t166 * t587;
t351 = g(1) * t571 + t168 * t587;
t347 = t292 * t374;
t346 = t294 * t373;
t345 = t296 * t372;
t344 = t286 * t347;
t343 = t290 * t345;
t342 = t288 * t346;
t335 = t275 * t46 + t374 * t621;
t334 = t276 * t46 + t374 * t622;
t333 = t275 * t47 + t373 * t621;
t332 = t276 * t47 + t373 * t622;
t331 = t275 * t48 + t372 * t621;
t330 = t276 * t48 + t372 * t622;
t301 = pkin(1) * g(1);
t214 = g(1) * qJ(2,1) + t586;
t213 = -g(2) * qJ(2,1) + t301;
t212 = g(1) * qJ(2,2) + t586;
t211 = -g(2) * qJ(2,2) + t301;
t210 = g(1) * qJ(2,3) + t586;
t209 = -g(2) * qJ(2,3) + t301;
t114 = t223 * t229 * t480 + t283 * t468;
t113 = t220 * t228 * t480 + t283 * t469;
t112 = t217 * t227 * t480 + t283 * t470;
t99 = t114 * t275 + t222 * t407;
t98 = t113 * t275 + t219 * t407;
t97 = t112 * t275 + t216 * t407;
t96 = -t114 * t276 + t222 * t408;
t95 = -t113 * t276 + t219 * t408;
t94 = -t112 * t276 + t216 * t408;
t78 = -t99 * t290 - t96 * t296;
t77 = -t98 * t288 - t95 * t294;
t76 = -t97 * t286 - t94 * t292;
t75 = -t96 * t290 + t99 * t296;
t74 = -t95 * t288 + t98 * t294;
t73 = -t94 * t286 + t97 * t292;
t72 = t114 * t594 + t372 * t623;
t71 = t113 * t595 + t373 * t623;
t70 = t112 * t596 + t374 * t623;
t42 = t369 + 0.2e1 * t532 - t108;
t41 = t370 + 0.2e1 * t531 - t107;
t40 = t371 + 0.2e1 * t530 - t106;
t39 = t48 * t594 + (t90 * t607 - t69 * t93) * t482;
t38 = t47 * t595 + (t89 * t607 - t68 * t92) * t484;
t37 = t46 * t596 + (t88 * t607 - t67 * t91) * t486;
t36 = t290 * t330 + t296 * t331;
t35 = t288 * t332 + t294 * t333;
t34 = t286 * t334 + t292 * t335;
t33 = t290 * t331 - t296 * t330;
t32 = t288 * t333 - t294 * t332;
t31 = t286 * t335 - t292 * t334;
t24 = g(1) * t135 + g(2) * t360 + 0.2e1 * t45 + t619;
t23 = g(1) * t134 + g(2) * t361 + 0.2e1 * t44 + t618;
t22 = g(1) * t133 + g(2) * t362 + 0.2e1 * t43 + t617;
t18 = -t45 - t629 * t254 + (-t423 + (-t90 * qJ(2,1) - t542 / 0.4e1 + t535 / 0.8e1) * t255) * t222 - t105;
t17 = -t44 - t628 * t251 + (-t424 + (-t89 * qJ(2,2) - t543 / 0.4e1 + t533 / 0.8e1) * t252) * t219 - t104;
t16 = -t43 - t627 * t248 + (-t425 + (-t88 * qJ(2,3) - t544 / 0.4e1 + t534 / 0.8e1) * t249) * t216 - t103;
t15 = (-0.4e1 * t342 + t47 - 0.2e1 * t516) * t264 + (t420 / 0.2e1 - t450 * t373) * t444 + 0.2e1 * t342 + t516;
t14 = (-0.4e1 * t343 + t48 - 0.2e1 * t515) * t264 + (t421 / 0.2e1 - t449 * t372) * t444 + 0.2e1 * t343 + t515;
t13 = (-0.4e1 * t344 + t46 - 0.2e1 * t517) * t264 + (t422 / 0.2e1 - t451 * t374) * t444 + 0.2e1 * t344 + t517;
t12 = (-t515 / 0.2e1 + t117 / 0.4e1 - t120 / 0.4e1 - t129 / 0.4e1 - t141 / 0.4e1) * t445 + t563 * t421 + ((-t57 * t255 / 0.8e1 + t256 * t397) * t445 - ((0.8e1 * t273 - 0.4e1) * t264 - 0.8e1 * t385 - 0.4e1 * t273 + 0.2e1) * t409) * t521;
t11 = (-t516 / 0.2e1 + t116 / 0.4e1 - t119 / 0.4e1 - t128 / 0.4e1 - t140 / 0.4e1) * t445 + t563 * t420 + ((-t56 * t252 / 0.8e1 + t253 * t398) * t445 - ((0.8e1 * t272 - 0.4e1) * t264 - 0.8e1 * t384 - 0.4e1 * t272 + 0.2e1) * t410) * t523;
t10 = (-t517 / 0.2e1 + t115 / 0.4e1 - t118 / 0.4e1 - t127 / 0.4e1 - t139 / 0.4e1) * t445 + t563 * t422 + ((-t55 * t249 / 0.8e1 + t250 * t399) * t445 - ((0.8e1 * t271 - 0.4e1) * t264 - 0.8e1 * t386 - 0.4e1 * t271 + 0.2e1) * t406) * t525;
t9 = (t291 * t213 - t214 * t297) * t246 + (t213 * t297 + t291 * t214) * t243 + (t369 + t532) * qJ(2,1) + (t45 + t619) * pkin(1);
t8 = (t289 * t211 - t212 * t295) * t245 + (t211 * t295 + t289 * t212) * t242 + (t370 + t531) * qJ(2,2) + (t44 + t618) * pkin(1);
t7 = (t287 * t209 - t210 * t293) * t244 + (t209 * t293 + t287 * t210) * t241 + (t371 + t530) * qJ(2,3) + (t43 + t617) * pkin(1);
t6 = t232 * t357 - t624 * t568 + t72 * t565 - t351 + t465 + (-t459 * t48 + t345) * pkin(2);
t5 = t231 * t358 - t625 * t569 + t71 * t566 - t352 + t466 + (-t460 * t47 + t346) * pkin(2);
t4 = t230 * t359 - t626 * t570 + t70 * t567 - t353 + t467 + (-t461 * t46 + t347) * pkin(2);
t3 = t48 * t560 + t226 * t357 + t624 * t565 + t72 * t568 - pkin(2) * (-t290 * t372 - t48 * t296) - t354 + t462;
t2 = t47 * t561 + t225 * t358 + t625 * t566 + t71 * t569 - pkin(2) * (-t288 * t373 - t47 * t294) - t355 + t463;
t1 = t46 * t562 + t224 * t359 + t626 * t567 + t70 * t570 - pkin(2) * (-t286 * t374 - t46 * t292) - t356 + t464;
t25 = [t46 * t505 + t47 * t504 + t48 * t503, t103 * t505 + t104 * t504 + t105 * t503, t106 * t505 + t107 * t504 + t108 * t503, (-(-t360 * t604 - t537) * t254 - (-t361 * t605 - t539) * t251 - (-t362 * t606 - t541) * t248) * t276, (-(t24 * t360 + t537) * t254 - (t23 * t361 + t539) * t251 - (t22 * t362 + t541) * t248) * t275, t40 * t505 + t41 * t504 + t42 * t503 + t428 * t84 + t433 * t83 + t438 * t82, -(t18 * t84 - t360 * t9) * t254 - (t17 * t83 - t361 * t8) * t251 - (t16 * t82 - t362 * t7) * t248, t13 * t505 + t14 * t503 + t15 * t504, t10 * t505 + t11 * t504 + t12 * t503, t503 * t75 + t504 * t74 + t505 * t73, t503 * t78 + t504 * t77 + t505 * t76, 0, -(-t3 * t360 + t33 * t84) * t254 - (-t2 * t361 + t32 * t83) * t251 - (-t1 * t362 + t31 * t82) * t248, -(t36 * t84 - t360 * t6) * t254 - (t35 * t83 - t361 * t5) * t251 - (t34 * t82 - t362 * t4) * t248, t285 - g(1); -t46 * t502 - t47 * t501 - t48 * t500, -t103 * t502 - t104 * t501 - t105 * t500, -t106 * t502 - t107 * t501 - t108 * t500, (-(t135 * t604 - t536) * t254 - (t134 * t605 - t538) * t251 - (t133 * t606 - t540) * t248) * t276, (-(-t135 * t24 + t536) * t254 - (-t134 * t23 + t538) * t251 - (-t133 * t22 + t540) * t248) * t275, -t40 * t502 - t41 * t501 - t42 * t500 + t428 * t87 + t433 * t86 + t438 * t85, -(t135 * t9 + t18 * t87) * t254 - (t134 * t8 + t17 * t86) * t251 - (t133 * t7 + t16 * t85) * t248, -t13 * t502 - t14 * t500 - t15 * t501, -t10 * t502 - t11 * t501 - t12 * t500, -t500 * t75 - t501 * t74 - t502 * t73, -t500 * t78 - t501 * t77 - t502 * t76, 0, -(t135 * t3 + t33 * t87) * t254 - (t134 * t2 + t32 * t86) * t251 - (t1 * t133 + t31 * t85) * t248, -(t135 * t6 + t36 * t87) * t254 - (t134 * t5 + t35 * t86) * t251 - (t133 * t4 + t34 * t85) * t248, t284 - g(2); -t412 * t48 - t414 * t47 - t416 * t46, -t103 * t416 - t104 * t414 - t105 * t412, -t106 * t416 - t107 * t414 - t108 * t412, (-(t21 * t598 - t527 / 0.2e1) * t483 - (t20 * t599 - t528 / 0.2e1) * t485 - (t19 * t600 - t529 / 0.2e1) * t487) * t276, (-(-t229 * t24 + t527 / 0.2e1) * t483 - (-t228 * t23 + t528 / 0.2e1) * t485 - (-t227 * t22 + t529 / 0.2e1) * t487) * t275, t217 * t520 * t579 + t220 * t519 * t578 + t223 * t518 * t577 - t40 * t416 - t41 * t414 - t412 * t42, -(t18 * t577 + t229 * t9) * t483 - (t17 * t578 + t228 * t8) * t485 - (t16 * t579 + t227 * t7) * t487, -t13 * t416 - t14 * t412 - t15 * t414 - 0.2e1 * ((t449 * t479 + t452 * t475) * t426 + (t450 * t479 + t452 * t473) * t431 + (t451 * t479 + t452 * t477) * t436) * t312, -t10 * t416 - t11 * t414 - t12 * t412 + ((-t273 * t563 + 0.4e1 * t385 + t564) * t426 + (-t272 * t563 + 0.4e1 * t384 + t564) * t431 + (-t271 * t563 + 0.4e1 * t386 + t564) * t436) * t312, -t73 * t416 - t74 * t414 - t75 * t412 + ((t296 * t275 + t290 * t276) * t522 + (t294 * t275 + t288 * t276) * t524 + (t292 * t275 + t286 * t276) * t526) * t312, -t76 * t416 - t77 * t414 - t78 * t412 + ((t472 - t476) * t522 + (-t288 * t275 + t294 * t276) * t524 + (t474 - t478) * t526) * t312, (t112 * t215 + t113 * t218 + t114 * t221) * t312, -t3 * t412 + t33 * t379 + (t39 * t568 - g(3) * t235 + (t226 / 0.2e1 + t290 / 0.2e1) * t393 + t354 + t462) * t468 - t2 * t414 + t32 * t381 + (t38 * t569 - g(3) * t234 + (t225 / 0.2e1 + t288 / 0.2e1) * t394 + t355 + t463) * t469 - t1 * t416 + t31 * t383 + (t37 * t570 - g(3) * t233 + (t224 / 0.2e1 + t286 / 0.2e1) * t395 + t356 + t464) * t470, -t6 * t412 + t36 * t379 + (t39 * t565 + g(3) * t229 + (t232 / 0.2e1 + t296 / 0.2e1) * t393 + t351 + t465) * t468 - t5 * t414 + t35 * t381 + (t38 * t566 + g(3) * t228 + (t231 / 0.2e1 + t294 / 0.2e1) * t394 + t352 + t466) * t469 - t4 * t416 + t34 * t383 + (t37 * t567 + g(3) * t227 + (t230 / 0.2e1 + t292 / 0.2e1) * t395 + t353 + t467) * t470, t283 - g(3);];
tauX_reg  = t25;
