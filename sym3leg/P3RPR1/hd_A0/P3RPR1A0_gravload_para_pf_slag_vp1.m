% Calculate Gravitation load for parallel robot
% P3RPR1A0
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
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPR1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:05
% EndTime: 2019-05-03 14:58:05
% DurationCPUTime: 0.26s
% Computational Cost: add. (309->86), mult. (492->161), div. (24->3), fcn. (362->14), ass. (0->77)
t480 = legFrame(3,3);
t469 = sin(t480);
t472 = cos(t480);
t457 = -g(1) * t469 + g(2) * t472;
t460 = g(1) * t472 + g(2) * t469;
t477 = rSges(2,3) + qJ(2,3);
t483 = sin(qJ(1,3));
t486 = cos(qJ(1,3));
t489 = pkin(1) + rSges(2,1);
t433 = ((-t457 * t489 - t460 * t477) * m(2) + m(1) * (-rSges(1,1) * t457 + rSges(1,2) * t460)) * t486 + t483 * ((-t457 * t477 + t460 * t489) * m(2) + m(1) * (rSges(1,1) * t460 + rSges(1,2) * t457));
t494 = 0.1e1 / qJ(2,3);
t508 = t433 * t494;
t481 = legFrame(2,3);
t470 = sin(t481);
t473 = cos(t481);
t458 = -g(1) * t470 + g(2) * t473;
t461 = g(1) * t473 + g(2) * t470;
t478 = rSges(2,3) + qJ(2,2);
t484 = sin(qJ(1,2));
t487 = cos(qJ(1,2));
t434 = ((-t458 * t489 - t461 * t478) * m(2) + m(1) * (-rSges(1,1) * t458 + rSges(1,2) * t461)) * t487 + t484 * ((-t458 * t478 + t461 * t489) * m(2) + m(1) * (rSges(1,1) * t461 + rSges(1,2) * t458));
t495 = 0.1e1 / qJ(2,2);
t507 = t434 * t495;
t482 = legFrame(1,3);
t471 = sin(t482);
t474 = cos(t482);
t459 = -g(1) * t471 + g(2) * t474;
t462 = g(1) * t474 + g(2) * t471;
t479 = rSges(2,3) + qJ(2,1);
t485 = sin(qJ(1,1));
t488 = cos(qJ(1,1));
t435 = ((-t459 * t489 - t462 * t479) * m(2) + m(1) * (-rSges(1,1) * t459 + rSges(1,2) * t462)) * t488 + t485 * ((-t459 * t479 + t462 * t489) * m(2) + m(1) * (rSges(1,1) * t462 + rSges(1,2) * t459));
t496 = 0.1e1 / qJ(2,1);
t506 = t435 * t496;
t442 = -t457 * t486 + t460 * t483;
t505 = t442 * t494;
t443 = -t458 * t487 + t461 * t484;
t504 = t443 * t495;
t444 = -t459 * t488 + t462 * t485;
t503 = t444 * t496;
t502 = koppelP(1,1);
t501 = koppelP(2,1);
t500 = koppelP(3,1);
t499 = koppelP(1,2);
t498 = koppelP(2,2);
t497 = koppelP(3,2);
t493 = rSges(3,1);
t492 = rSges(3,2);
t491 = xP(3);
t490 = pkin(1) + pkin(2);
t476 = cos(t491);
t475 = sin(t491);
t468 = qJ(2,1) * t485 + t488 * t490;
t467 = qJ(2,2) * t484 + t487 * t490;
t466 = qJ(2,3) * t483 + t486 * t490;
t465 = -qJ(2,1) * t488 + t485 * t490;
t464 = -qJ(2,2) * t487 + t484 * t490;
t463 = -qJ(2,3) * t486 + t483 * t490;
t456 = -t475 * t499 + t476 * t502;
t455 = -t475 * t498 + t476 * t501;
t454 = -t475 * t497 + t476 * t500;
t453 = -t475 * t502 - t476 * t499;
t452 = -t475 * t501 - t476 * t498;
t451 = -t475 * t500 - t476 * t497;
t450 = -t471 * t485 + t474 * t488;
t449 = t471 * t488 + t474 * t485;
t448 = -t470 * t484 + t473 * t487;
t447 = t470 * t487 + t473 * t484;
t446 = -t469 * t483 + t472 * t486;
t445 = t469 * t486 + t472 * t483;
t441 = -t465 * t471 + t468 * t474;
t440 = -t464 * t470 + t467 * t473;
t439 = -t463 * t469 + t466 * t472;
t438 = t465 * t474 + t468 * t471;
t437 = t464 * t473 + t467 * t470;
t436 = t463 * t472 + t466 * t469;
t1 = [t446 * t508 + t448 * t507 + t450 * t506 - m(3) * g(1) + (-t439 * t505 - t440 * t504 - t441 * t503) * m(2); t445 * t508 + t447 * t507 + t449 * t506 - m(3) * g(2) + (-t436 * t505 - t437 * t504 - t438 * t503) * m(2); m(3) * ((g(1) * t493 + g(2) * t492) * t475 + (g(1) * t492 - g(2) * t493) * t476) + ((t449 * t456 + t450 * t453) * t435 - (t438 * t456 + t441 * t453) * m(2) * t444) * t496 + ((t447 * t455 + t448 * t452) * t434 - (t437 * t455 + t440 * t452) * m(2) * t443) * t495 + ((t445 * t454 + t446 * t451) * t433 - (t436 * t454 + t439 * t451) * m(2) * t442) * t494;];
taugX  = t1;
