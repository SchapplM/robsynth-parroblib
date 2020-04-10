% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR2G2P3A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G2P3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR2G2P3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:47
% EndTime: 2020-03-09 21:21:47
% DurationCPUTime: 0.48s
% Computational Cost: add. (318->93), mult. (449->170), div. (126->5), fcn. (492->24), ass. (0->87)
t1488 = qJ(2,3) + qJ(3,3);
t1473 = sin(t1488);
t1494 = sin(qJ(3,3));
t1485 = 0.1e1 / t1494;
t1495 = sin(qJ(2,3));
t1534 = (t1495 * pkin(1) + pkin(2) * t1473) * t1485;
t1489 = qJ(2,2) + qJ(3,2);
t1474 = sin(t1489);
t1496 = sin(qJ(3,2));
t1486 = 0.1e1 / t1496;
t1497 = sin(qJ(2,2));
t1533 = (t1497 * pkin(1) + pkin(2) * t1474) * t1486;
t1490 = qJ(2,1) + qJ(3,1);
t1475 = sin(t1490);
t1498 = sin(qJ(3,1));
t1487 = 0.1e1 / t1498;
t1499 = sin(qJ(2,1));
t1532 = (t1499 * pkin(1) + pkin(2) * t1475) * t1487;
t1531 = t1473 * t1485;
t1530 = t1474 * t1486;
t1529 = t1475 * t1487;
t1491 = legFrame(3,2);
t1479 = sin(t1491);
t1528 = t1479 * t1485;
t1492 = legFrame(2,2);
t1480 = sin(t1492);
t1527 = t1480 * t1486;
t1493 = legFrame(1,2);
t1481 = sin(t1493);
t1526 = t1481 * t1487;
t1482 = cos(t1491);
t1525 = t1482 * t1485;
t1483 = cos(t1492);
t1524 = t1483 * t1486;
t1484 = cos(t1493);
t1523 = t1484 * t1487;
t1522 = t1495 * t1494;
t1521 = t1497 * t1496;
t1520 = t1499 * t1498;
t1500 = cos(qJ(3,3));
t1501 = cos(qJ(2,3));
t1458 = (pkin(2) * t1500 + pkin(1)) * t1501 - pkin(2) * t1522;
t1519 = t1458 * t1528;
t1518 = t1458 * t1525;
t1502 = cos(qJ(3,2));
t1503 = cos(qJ(2,2));
t1459 = (pkin(2) * t1502 + pkin(1)) * t1503 - pkin(2) * t1521;
t1517 = t1459 * t1527;
t1516 = t1459 * t1524;
t1504 = cos(qJ(3,1));
t1505 = cos(qJ(2,1));
t1460 = (pkin(2) * t1504 + pkin(1)) * t1505 - pkin(2) * t1520;
t1515 = t1460 * t1526;
t1514 = t1460 * t1523;
t1461 = t1501 * t1500 - t1522;
t1513 = t1461 * t1528;
t1512 = t1461 * t1525;
t1462 = t1503 * t1502 - t1521;
t1511 = t1462 * t1527;
t1510 = t1462 * t1524;
t1463 = t1505 * t1504 - t1520;
t1509 = t1463 * t1526;
t1508 = t1463 * t1523;
t1507 = 0.1e1 / pkin(1);
t1506 = 0.1e1 / pkin(2);
t1478 = cos(t1490);
t1477 = cos(t1489);
t1476 = cos(t1488);
t1469 = t1484 * g(1) - t1481 * g(2);
t1468 = t1483 * g(1) - t1480 * g(2);
t1467 = t1482 * g(1) - t1479 * g(2);
t1466 = t1481 * g(1) + t1484 * g(2);
t1465 = t1480 * g(1) + t1483 * g(2);
t1464 = t1479 * g(1) + t1482 * g(2);
t1457 = g(3) * t1505 + t1466 * t1499;
t1456 = g(3) * t1503 + t1465 * t1497;
t1455 = g(3) * t1501 + t1464 * t1495;
t1454 = -g(3) * t1499 + t1466 * t1505;
t1453 = -g(3) * t1497 + t1465 * t1503;
t1452 = -g(3) * t1495 + t1464 * t1501;
t1451 = g(3) * t1478 + t1466 * t1475;
t1450 = g(3) * t1477 + t1465 * t1474;
t1449 = g(3) * t1476 + t1464 * t1473;
t1448 = -g(3) * t1475 + t1466 * t1478;
t1447 = -g(3) * t1474 + t1465 * t1477;
t1446 = -g(3) * t1473 + t1464 * t1476;
t1 = [(-t1482 * t1467 - t1483 * t1468 - t1484 * t1469) * MDP(1) - g(1) * MDP(8) + ((t1455 * t1513 + t1456 * t1511 + t1457 * t1509) * MDP(3) + (t1452 * t1513 + t1453 * t1511 + t1454 * t1509) * MDP(4) + (t1449 * t1513 + t1450 * t1511 + t1451 * t1509) * MDP(6) + (t1446 * t1513 + t1447 * t1511 + t1448 * t1509) * MDP(7) + ((-t1449 * t1519 - t1450 * t1517 - t1451 * t1515) * MDP(6) + (-t1446 * t1519 - t1447 * t1517 - t1448 * t1515) * MDP(7)) * t1506) * t1507; (t1479 * t1467 + t1480 * t1468 + t1481 * t1469) * MDP(1) - g(2) * MDP(8) + ((t1455 * t1512 + t1456 * t1510 + t1457 * t1508) * MDP(3) + (t1452 * t1512 + t1453 * t1510 + t1454 * t1508) * MDP(4) + (t1449 * t1512 + t1450 * t1510 + t1451 * t1508) * MDP(6) + (t1446 * t1512 + t1447 * t1510 + t1448 * t1508) * MDP(7) + ((-t1449 * t1518 - t1450 * t1516 - t1451 * t1514) * MDP(6) + (-t1446 * t1518 - t1447 * t1516 - t1448 * t1514) * MDP(7)) * t1506) * t1507; -g(3) * MDP(8) + ((-t1455 * t1531 - t1456 * t1530 - t1457 * t1529) * MDP(3) + (-t1452 * t1531 - t1453 * t1530 - t1454 * t1529) * MDP(4) + (-t1449 * t1531 - t1450 * t1530 - t1451 * t1529) * MDP(6) + (-t1446 * t1531 - t1447 * t1530 - t1448 * t1529) * MDP(7) + ((t1449 * t1534 + t1450 * t1533 + t1451 * t1532) * MDP(6) + (t1446 * t1534 + t1447 * t1533 + t1448 * t1532) * MDP(7)) * t1506) * t1507;];
taugX  = t1;
