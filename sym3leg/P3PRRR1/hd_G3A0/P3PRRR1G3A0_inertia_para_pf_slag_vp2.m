% Calculate inertia matrix for parallel robot
% P3PRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:50
% EndTime: 2020-03-09 21:06:50
% DurationCPUTime: 0.53s
% Computational Cost: add. (1974->97), mult. (1212->195), div. (405->5), fcn. (1179->24), ass. (0->112)
t306 = 0.1e1 / pkin(2);
t293 = pkin(7) + qJ(2,1);
t284 = qJ(3,1) + t293;
t272 = sin(t284);
t275 = cos(t284);
t278 = sin(t293);
t281 = cos(t293);
t358 = 0.1e1 / (t272 * t281 - t278 * t275);
t345 = t358 * t306;
t292 = pkin(7) + qJ(2,2);
t283 = qJ(3,2) + t292;
t271 = sin(t283);
t274 = cos(t283);
t277 = sin(t292);
t280 = cos(t292);
t359 = 0.1e1 / (t271 * t280 - t277 * t274);
t349 = t359 * t306;
t291 = pkin(7) + qJ(2,3);
t282 = qJ(3,3) + t291;
t270 = sin(t282);
t273 = cos(t282);
t276 = sin(t291);
t279 = cos(t291);
t360 = 0.1e1 / (t270 * t279 - t276 * t273);
t353 = t360 * t306;
t297 = legFrame(1,2);
t290 = cos(t297);
t372 = t290 * t345;
t296 = legFrame(2,2);
t289 = cos(t296);
t371 = t289 * t349;
t295 = legFrame(3,2);
t288 = cos(t295);
t370 = t288 * t353;
t307 = (mrSges(3,1) * cos(qJ(3,1)) - mrSges(3,2) * sin(qJ(3,1))) * pkin(2);
t332 = m(3) * pkin(2) ^ 2 + Ifges(2,3) + Ifges(3,3);
t257 = 0.2e1 * t307 + t332;
t263 = pkin(2) * t281 + pkin(3) * t275;
t269 = Ifges(3,3) + t307;
t305 = 0.1e1 / pkin(3);
t336 = t269 * t305;
t369 = -t257 * t275 + t263 * t336;
t308 = (mrSges(3,1) * cos(qJ(3,2)) - mrSges(3,2) * sin(qJ(3,2))) * pkin(2);
t256 = 0.2e1 * t308 + t332;
t262 = pkin(2) * t280 + pkin(3) * t274;
t268 = Ifges(3,3) + t308;
t338 = t268 * t305;
t368 = -t256 * t274 + t262 * t338;
t309 = (mrSges(3,1) * cos(qJ(3,3)) - mrSges(3,2) * sin(qJ(3,3))) * pkin(2);
t255 = 0.2e1 * t309 + t332;
t261 = pkin(2) * t279 + pkin(3) * t273;
t267 = Ifges(3,3) + t309;
t340 = t267 * t305;
t367 = -t255 * t273 + t261 * t340;
t357 = Ifges(3,3) * t305;
t366 = t263 * t357 - t269 * t275;
t365 = t262 * t357 - t268 * t274;
t364 = t261 * t357 - t267 * t273;
t363 = t288 * t360;
t362 = t289 * t359;
t361 = t290 * t358;
t258 = pkin(2) * t276 + pkin(3) * t270;
t356 = t360 * t258;
t355 = t360 * t270;
t285 = sin(t295);
t354 = t360 * t285;
t259 = pkin(2) * t277 + pkin(3) * t271;
t352 = t359 * t259;
t351 = t359 * t271;
t286 = sin(t296);
t350 = t359 * t286;
t260 = pkin(2) * t278 + pkin(3) * t272;
t348 = t358 * t260;
t347 = t358 * t272;
t287 = sin(t297);
t346 = t358 * t287;
t328 = t261 * t354;
t327 = t273 * t354;
t326 = t273 * t363;
t325 = t285 * t353;
t324 = t262 * t350;
t323 = t274 * t350;
t322 = t274 * t362;
t321 = t286 * t349;
t320 = t263 * t346;
t319 = t275 * t346;
t318 = t275 * t361;
t317 = t287 * t345;
t316 = t261 * t363;
t315 = t262 * t362;
t314 = t263 * t361;
t294 = m(1) + m(2) + m(3);
t310 = (t285 * t288 + t286 * t289 + t287 * t290) * t294;
t245 = (t260 * t357 - t269 * t272) * t345;
t244 = (t259 * t357 - t268 * t271) * t349;
t243 = (t258 * t357 - t267 * t270) * t353;
t242 = t366 * t372;
t241 = t365 * t371;
t240 = t364 * t370;
t239 = t366 * t317;
t238 = t365 * t321;
t237 = t364 * t325;
t236 = (-t257 * t272 + t260 * t336) * t345;
t235 = (-t256 * t271 + t259 * t338) * t349;
t234 = (-t255 * t270 + t258 * t340) * t353;
t233 = t369 * t372;
t232 = t368 * t371;
t231 = t367 * t370;
t230 = t369 * t317;
t229 = t368 * t321;
t228 = t367 * t325;
t1 = [m(4) + (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) * t294 + (-t231 * t326 - t232 * t322 - t233 * t318 + (t240 * t316 + t241 * t315 + t242 * t314) * t305) * t306, (t231 * t327 + t232 * t323 + t233 * t319 + (-t240 * t328 - t241 * t324 - t242 * t320) * t305) * t306 + t310, (t231 * t355 + t232 * t351 + t233 * t347 + (-t240 * t356 - t241 * t352 - t242 * t348) * t305) * t306; (t228 * t326 + t229 * t322 + t230 * t318 + (-t237 * t316 - t238 * t315 - t239 * t314) * t305) * t306 + t310, m(4) + (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) * t294 + (-t228 * t327 - t229 * t323 - t230 * t319 + (t237 * t328 + t238 * t324 + t239 * t320) * t305) * t306, (-t228 * t355 - t229 * t351 - t230 * t347 + (t237 * t356 + t238 * t352 + t239 * t348) * t305) * t306; (t234 * t326 + t235 * t322 + t236 * t318 + (-t243 * t316 - t244 * t315 - t245 * t314) * t305) * t306, (-t234 * t327 - t235 * t323 - t236 * t319 + (t243 * t328 + t244 * t324 + t245 * t320) * t305) * t306, m(4) + (-t234 * t355 - t235 * t351 - t236 * t347 + (t243 * t356 + t244 * t352 + t245 * t348) * t305) * t306;];
MX  = t1;
