package ru.processtech.forome.annotation.config.connector;

import net.minidev.json.JSONObject;
import ru.processtech.forome.annotation.config.connector.database.DatabaseConfigConnector;

public class GTFConfigConfigConnector extends DatabaseConfigConnector {

	public GTFConfigConfigConnector(JSONObject parse) {
		super(parse);
	}
}
