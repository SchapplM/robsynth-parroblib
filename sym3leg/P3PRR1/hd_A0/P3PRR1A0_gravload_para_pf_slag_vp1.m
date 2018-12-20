% Calculate Gravitation load for parallel robot
% P3PRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m [3x1]
%   mass of all robot links (including platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:42
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3PRR1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:42:35
% EndTime: 2018-12-20 17:42:35
% DurationCPUTime: 0.19s
% Computational Cost: add. (132->55), mult. (249->105), div. (36->4), fcn. (218->14), ass. (0->56)
t402 = m(2) / pkin(2);
t376 = legFrame(3,3);
t365 = sin(t376);
t368 = cos(t376);
t359 = -g(1) * t365 + g(2) * t368;
t362 = g(1) * t368 + g(2) * t365;
t379 = sin(qJ(2,3));
t373 = 0.1e1 / t379;
t382 = cos(qJ(2,3));
t401 = (t379 * (rSges(2,1) * t362 + rSges(2,2) * t359) - t382 * (rSges(2,1) * t359 - rSges(2,2) * t362)) * t373;
t377 = legFrame(2,3);
t366 = sin(t377);
t369 = cos(t377);
t360 = -g(1) * t366 + g(2) * t369;
t363 = g(1) * t369 + g(2) * t366;
t380 = sin(qJ(2,2));
t374 = 0.1e1 / t380;
t383 = cos(qJ(2,2));
t400 = (t380 * (rSges(2,1) * t363 + rSges(2,2) * t360) - t383 * (rSges(2,1) * t360 - rSges(2,2) * t363)) * t374;
t378 = legFrame(1,3);
t367 = sin(t378);
t370 = cos(t378);
t361 = -g(1) * t367 + g(2) * t370;
t364 = g(1) * t370 + g(2) * t367;
t381 = sin(qJ(2,1));
t375 = 0.1e1 / t381;
t384 = cos(qJ(2,1));
t399 = (t381 * (rSges(2,1) * t364 + rSges(2,2) * t361) - t384 * (rSges(2,1) * t361 - rSges(2,2) * t364)) * t375;
t398 = t359 * t373;
t397 = t360 * t374;
t396 = t361 * t375;
t394 = koppelP(1,1);
t393 = koppelP(2,1);
t392 = koppelP(3,1);
t391 = koppelP(1,2);
t390 = koppelP(2,2);
t389 = koppelP(3,2);
t388 = rSges(3,1);
t387 = rSges(3,2);
t386 = xP(3);
t385 = m(1) + m(2);
t372 = cos(t386);
t371 = sin(t386);
t358 = -t371 * t391 + t372 * t394;
t357 = -t371 * t390 + t372 * t393;
t356 = -t371 * t389 + t372 * t392;
t355 = -t371 * t394 - t372 * t391;
t354 = -t371 * t393 - t372 * t390;
t353 = -t371 * t392 - t372 * t389;
t352 = t367 * t384 + t370 * t381;
t351 = -t367 * t381 + t370 * t384;
t350 = t366 * t383 + t369 * t380;
t349 = -t366 * t380 + t369 * t383;
t348 = t365 * t382 + t368 * t379;
t347 = -t365 * t379 + t368 * t382;
t1 = [-m(3) * g(1) + (-t347 * t398 - t349 * t397 - t351 * t396) * t385 + (-t368 * t401 - t369 * t400 - t370 * t399) * t402; -m(3) * g(2) + (-t348 * t398 - t350 * t397 - t352 * t396) * t385 + (-t365 * t401 - t366 * t400 - t367 * t399) * t402; m(3) * ((g(1) * t388 + g(2) * t387) * t371 + (g(1) * t387 - g(2) * t388) * t372) + (-(t351 * t355 + t352 * t358) * t396 - (t349 * t354 + t350 * t357) * t397 - (t347 * t353 + t348 * t356) * t398) * t385 + ((-t355 * t370 - t358 * t367) * t399 + (-t354 * t369 - t357 * t366) * t400 + (-t353 * t368 - t356 * t365) * t401) * t402;];
taugX  = t1;
