% Calculate inertia matrix for parallel robot
% P3PRRR1G2A0
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
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:22
% EndTime: 2020-03-09 21:18:22
% DurationCPUTime: 0.59s
% Computational Cost: add. (1974->97), mult. (1212->194), div. (405->5), fcn. (1179->24), ass. (0->108)
t308 = legFrame(3,2);
t301 = cos(t308);
t304 = pkin(7) + qJ(2,3);
t295 = qJ(3,3) + t304;
t283 = sin(t295);
t286 = cos(t295);
t289 = sin(t304);
t292 = cos(t304);
t372 = 0.1e1 / (-t283 * t292 + t286 * t289);
t375 = t301 * t372;
t309 = legFrame(2,2);
t302 = cos(t309);
t305 = pkin(7) + qJ(2,2);
t296 = qJ(3,2) + t305;
t284 = sin(t296);
t287 = cos(t296);
t290 = sin(t305);
t293 = cos(t305);
t371 = 0.1e1 / (-t284 * t293 + t287 * t290);
t374 = t302 * t371;
t310 = legFrame(1,2);
t303 = cos(t310);
t306 = pkin(7) + qJ(2,1);
t297 = qJ(3,1) + t306;
t285 = sin(t297);
t288 = cos(t297);
t291 = sin(t306);
t294 = cos(t306);
t370 = 0.1e1 / (-t285 * t294 + t288 * t291);
t373 = t303 * t370;
t319 = 0.1e1 / pkin(2);
t384 = t319 * t373;
t383 = t319 * t374;
t382 = t319 * t375;
t320 = (mrSges(3,1) * cos(qJ(3,1)) - mrSges(3,2) * sin(qJ(3,1))) * pkin(2);
t348 = m(3) * pkin(2) ^ 2 + Ifges(2,3) + Ifges(3,3);
t270 = 0.2e1 * t320 + t348;
t282 = Ifges(3,3) + t320;
t273 = pkin(2) * t291 + pkin(3) * t285;
t318 = 0.1e1 / pkin(3);
t355 = t273 * t318;
t381 = -t270 * t285 + t282 * t355;
t321 = (mrSges(3,1) * cos(qJ(3,2)) - mrSges(3,2) * sin(qJ(3,2))) * pkin(2);
t269 = 0.2e1 * t321 + t348;
t281 = Ifges(3,3) + t321;
t272 = pkin(2) * t290 + pkin(3) * t284;
t356 = t272 * t318;
t380 = -t269 * t284 + t281 * t356;
t322 = (mrSges(3,1) * cos(qJ(3,3)) - mrSges(3,2) * sin(qJ(3,3))) * pkin(2);
t268 = 0.2e1 * t322 + t348;
t280 = Ifges(3,3) + t322;
t271 = pkin(2) * t289 + pkin(3) * t283;
t357 = t271 * t318;
t379 = -t268 * t283 + t280 * t357;
t378 = Ifges(3,3) * t355 - t282 * t285;
t377 = Ifges(3,3) * t356 - t281 * t284;
t376 = Ifges(3,3) * t357 - t280 * t283;
t369 = t372 * (pkin(2) * t292 + pkin(3) * t286);
t298 = sin(t308);
t368 = t372 * t298;
t367 = t371 * (pkin(2) * t293 + pkin(3) * t287);
t299 = sin(t309);
t366 = t371 * t299;
t365 = t370 * (pkin(2) * t294 + pkin(3) * t288);
t300 = sin(t310);
t364 = t370 * t300;
t363 = t372 * t286;
t362 = t371 * t287;
t361 = t370 * t288;
t344 = t271 * t368;
t343 = t271 * t375;
t342 = t318 * t369;
t341 = t283 * t368;
t340 = t319 * t368;
t339 = t272 * t366;
t338 = t272 * t374;
t337 = t318 * t367;
t336 = t284 * t366;
t335 = t319 * t366;
t334 = t273 * t364;
t333 = t273 * t373;
t332 = t318 * t365;
t331 = t285 * t364;
t330 = t319 * t364;
t329 = t283 * t375;
t328 = t284 * t374;
t327 = t285 * t373;
t307 = m(1) + m(2) + m(3);
t323 = (t298 * t301 + t299 * t302 + t300 * t303) * t307;
t258 = (Ifges(3,3) * t332 - t282 * t361) * t319;
t257 = (Ifges(3,3) * t337 - t281 * t362) * t319;
t256 = (Ifges(3,3) * t342 - t280 * t363) * t319;
t255 = t378 * t330;
t254 = t377 * t335;
t253 = t376 * t340;
t252 = t378 * t384;
t251 = t377 * t383;
t250 = t376 * t382;
t249 = (-t270 * t361 + t282 * t332) * t319;
t248 = (-t269 * t362 + t281 * t337) * t319;
t247 = (-t268 * t363 + t280 * t342) * t319;
t246 = t381 * t330;
t245 = t380 * t335;
t244 = t379 * t340;
t243 = t381 * t384;
t242 = t380 * t383;
t241 = t379 * t382;
t1 = [m(4) + (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) * t307 + (-t241 * t329 - t242 * t328 - t243 * t327 + (t250 * t343 + t251 * t338 + t252 * t333) * t318) * t319, (t241 * t341 + t242 * t336 + t243 * t331 + (-t250 * t344 - t251 * t339 - t252 * t334) * t318) * t319 + t323, (-t241 * t363 - t242 * t362 - t243 * t361 + (t250 * t369 + t251 * t367 + t252 * t365) * t318) * t319; (t244 * t329 + t245 * t328 + t246 * t327 + (-t253 * t343 - t254 * t338 - t255 * t333) * t318) * t319 + t323, m(4) + (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) * t307 + (-t244 * t341 - t245 * t336 - t246 * t331 + (t253 * t344 + t254 * t339 + t255 * t334) * t318) * t319, (t244 * t363 + t245 * t362 + t246 * t361 + (-t253 * t369 - t254 * t367 - t255 * t365) * t318) * t319; (-t247 * t329 - t248 * t328 - t249 * t327 + (t256 * t343 + t257 * t338 + t258 * t333) * t318) * t319, (t247 * t341 + t248 * t336 + t249 * t331 + (-t256 * t344 - t257 * t339 - t258 * t334) * t318) * t319, m(4) + (-t247 * t363 - t248 * t362 - t249 * t361 + (t256 * t369 + t257 * t367 + t258 * t365) * t318) * t319;];
MX  = t1;
