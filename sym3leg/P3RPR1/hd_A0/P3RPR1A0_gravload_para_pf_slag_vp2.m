% Calculate Gravitation load for parallel robot
% P3RPR1G1P1A0
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
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPR1G1P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:08
% EndTime: 2019-05-03 14:58:08
% DurationCPUTime: 0.24s
% Computational Cost: add. (273->76), mult. (383->140), div. (24->3), fcn. (290->14), ass. (0->78)
t483 = mrSges(2,3) - mrSges(1,2);
t455 = legFrame(3,3);
t447 = sin(t455);
t450 = cos(t455);
t431 = -t447 * g(1) + t450 * g(2);
t434 = t450 * g(1) + t447 * g(2);
t443 = m(2) * qJ(2,3) + t483;
t446 = m(2) * pkin(1) + mrSges(1,1) + mrSges(2,1);
t458 = sin(qJ(1,3));
t461 = cos(qJ(1,3));
t407 = (-t431 * t446 - t443 * t434) * t461 - (t431 * t443 - t446 * t434) * t458;
t468 = 0.1e1 / qJ(2,3);
t482 = t407 * t468;
t456 = legFrame(2,3);
t448 = sin(t456);
t451 = cos(t456);
t432 = -t448 * g(1) + t451 * g(2);
t435 = t451 * g(1) + t448 * g(2);
t444 = m(2) * qJ(2,2) + t483;
t459 = sin(qJ(1,2));
t462 = cos(qJ(1,2));
t408 = (-t432 * t446 - t444 * t435) * t462 - (t432 * t444 - t446 * t435) * t459;
t469 = 0.1e1 / qJ(2,2);
t481 = t408 * t469;
t457 = legFrame(1,3);
t449 = sin(t457);
t452 = cos(t457);
t433 = -t449 * g(1) + t452 * g(2);
t436 = t452 * g(1) + t449 * g(2);
t445 = m(2) * qJ(2,1) + t483;
t460 = sin(qJ(1,1));
t463 = cos(qJ(1,1));
t409 = (-t433 * t446 - t445 * t436) * t463 - (t433 * t445 - t446 * t436) * t460;
t470 = 0.1e1 / qJ(2,1);
t480 = t409 * t470;
t416 = -t461 * t431 + t458 * t434;
t479 = t416 * t468;
t417 = -t462 * t432 + t459 * t435;
t478 = t417 * t469;
t418 = -t463 * t433 + t460 * t436;
t477 = t418 * t470;
t476 = koppelP(1,1);
t475 = koppelP(2,1);
t474 = koppelP(3,1);
t473 = koppelP(1,2);
t472 = koppelP(2,2);
t471 = koppelP(3,2);
t467 = mrSges(3,1);
t466 = mrSges(3,2);
t465 = xP(3);
t464 = pkin(1) + pkin(2);
t454 = cos(t465);
t453 = sin(t465);
t442 = t460 * qJ(2,1) + t464 * t463;
t441 = t459 * qJ(2,2) + t464 * t462;
t440 = t458 * qJ(2,3) + t464 * t461;
t439 = -t463 * qJ(2,1) + t460 * t464;
t438 = -t462 * qJ(2,2) + t459 * t464;
t437 = -t461 * qJ(2,3) + t458 * t464;
t430 = -t453 * t473 + t454 * t476;
t429 = -t453 * t472 + t454 * t475;
t428 = -t453 * t471 + t454 * t474;
t427 = -t453 * t476 - t454 * t473;
t426 = -t453 * t475 - t454 * t472;
t425 = -t453 * t474 - t454 * t471;
t424 = -t449 * t460 + t452 * t463;
t423 = t449 * t463 + t452 * t460;
t422 = -t448 * t459 + t451 * t462;
t421 = t448 * t462 + t451 * t459;
t420 = -t447 * t458 + t450 * t461;
t419 = t447 * t461 + t450 * t458;
t415 = -t449 * t439 + t442 * t452;
t414 = -t448 * t438 + t441 * t451;
t413 = -t447 * t437 + t440 * t450;
t412 = t439 * t452 + t449 * t442;
t411 = t438 * t451 + t448 * t441;
t410 = t437 * t450 + t447 * t440;
t1 = [t420 * t482 + t422 * t481 + t424 * t480 - g(1) * m(3) + (-t413 * t479 - t414 * t478 - t415 * t477) * m(2); t419 * t482 + t421 * t481 + t423 * t480 - g(2) * m(3) + (-t410 * t479 - t411 * t478 - t412 * t477) * m(2); -(-g(1) * t467 - g(2) * t466) * t453 + t454 * (g(1) * t466 - g(2) * t467) + ((t423 * t430 + t424 * t427) * t409 - (t412 * t430 + t415 * t427) * m(2) * t418) * t470 + ((t421 * t429 + t422 * t426) * t408 - (t411 * t429 + t414 * t426) * m(2) * t417) * t469 + ((t419 * t428 + t420 * t425) * t407 - (t410 * t428 + t413 * t425) * m(2) * t416) * t468;];
taugX  = t1;
